#Unscented Kalman Filter Implimentation for Cornell Rocketry Real-Time State Estimation
#Developed by Sam Noles in 2022

#PyPlot Plotting Examples: https://gist.github.com/gizmaa/7214002#multiaxis

#########info##############

#system state: z = [r v q w]
#r: position of R (tip of nose cone) wrt to O (ground) in inertial frame 1,2,3
#v: inertial frame derivative of r in intertial frame coordinates 4,5,6
#q: quaternion representing rotation from inertial to body frame 7,8,9,10
    #aka apply q to a vector written in I cordinates to get the same vector written in body frame coordinates
#w:IwB rotation of B frame wrt to the inertial frame 11,12,13

#estimator state: ztilda = [r v ds w]
#r: position of R (tip of nose cone) wrt to O (ground) in inertial frame 1,2,3
#v: inertial frame derivative of r in intertial frame coordinates 4,5,6
#ds: generalize Rodrigues error vector representing error from the quaternion in the system state
#w:IwB rotation of B frame wrt to the inertial frame 11,12,13


###Estimator#####

#y: Measurement 11x1: [y_accel, y_gyro, y_barox3, y_magb]
#R: Sensor noise covariance 10x10
#w: [ww1, ww2, ww3, wt, wF1, wF2, wF3, wm1, wm2, wm3] 
#Q: Process noise covariance 4x4

include("6dof.jl")

const global STATE_SIZE = 12 #estimator state size
const global W_SIZE = 10;
const global Y_SIZE = 10;
const global Be_I = [0;1.0;0.0] #direction of the Earths magnetic field in the inertial frame


#functions 

#implimentation of Gw(t,z)
function getGw(t::Float64, z::Vector{Float64}, dt::Float64, simParam::sim)
    #t: current time
    #z: state vector
    #dt: time step
    #rocket: rocket object that has all the rocket info
    #return: Gw STATE_SIZE x 4 matrix w = [ww1, ww2, ww3, wt]

    Gw = zeros(STATE_SIZE, W_SIZE);
    m = sum(simParam.rocket.massData.currentMass)

    #contribution from thrust uncertainty
    Gw[6, 4] = dt / m
    Gw[4:6, 5:7] = (dt / m) * diagm(ones(3))  #contribution of general forcing uncertainty
    Gw[10:12, 8:10] = dt * inv(getIg(simParam.rocket.massData)) #contribution of general moment uncertainty

    vAOI_I = getWind(t, z[3], simParam.simInputs.windData)
    vRAI_I = getVRA(z[4:6], vAOI_I)
    mag_vRAI_I = norm(vRAI_I)
    aoa = calcAoA(vRAI_I, z[7:10])
    mach = calcMach(vRAI_I, z[3])


    #wind uncertainty
    dfd_dwd = dt * 0.5 * getCd(aoa, mach, simParam.rocket.aeroData) * getA(aoa, simParam.rocket.aeroData) * expAtm(z[3]) * (vRAI_I * transpose(vRAI_I) / mag_vRAI_I + (I * mag_vRAI_I))/m

    Gw[4:6, 1:3] = dfd_dwd;
    
    return Gw

end

function getQ(t::Float64, z::Vector{Float64}, lv::rocket)

    Qwind = getWind(t, z[3])[2]
    thrust = motorThrustMass(t, lv.motorData, lv.massData.initalMass[2])[1]
    Ig_B = getIg(lv.massData)
    thrustVariation = getThrustVar(t)

    addativeForceProcessNoise = 625.0 #Newtons^2
    addativeMomentProcessNoise = 100.0  #(Nm)^2

    Q = zeros(W_SIZE,W_SIZE)
    Q[1:3,1:3] = Qwind
    Q[4,4] = thrustVariation * thrust
    Q[5:7,5:7] = addativeForceProcessNoise * diagm(ones(3))
    Q[8:10,8:10] = addativeMomentProcessNoise * diagm(ones(3))

    return Q
end

function sigmaPoints_nsigma(zhat::Vector{Float64}, n::Int, nsigma::Float64, Pkk::Matrix{Float64})
    #zhat: state estimate
    #n: size of state vector
    #nsigma: number of standard deviations to place the sigma point from the mean
    #Pkk: covariance of state estimate
    #returns: chi (n x n*2+1) matrix with columns as sigma points

    chi = zeros(n, n * 2 + 1)

    Phalf = cholesky(Pkk);

    chi[:, 1] = zhat #first sigma point is just the mean
    for i = 1:n
        
        chi[:, 1 + i] = zhat + nsigma * Phalf.L[:,i]
        chi[:, 1 + n + i] = zhat - nsigma * Phalf.L[:,i]

    end

    return chi
    
end

function nsigma_weights(nsigma::Float64)
    #nsigma: number of std the sigma points are spaced from the mean
    #returns: the three different weights to be used W0S, W0C, wi

    W0S = (nsigma^2 - STATE_SIZE)/nsigma^2
    W0C = W0S - nsigma^2/STATE_SIZE  + 3
    wi = 0.5 * 1/nsigma^2

    return W0S, W0C, wi


end

function sigmaPoints(zhat::Vector{Float64}, n::Int, Wo::Float64, Pkk::Matrix{Float64})
    #zhat: state estimate
    #n: size of state vector
    #Wo: weighting of first sigma point
    #Pkk: covariance of state estimate
    #returns: chi (n x n*2+1) matrix with columns as sigma points

    nsigma = sqrt(n/(1-Wo))

    return sigmaPoints_nsigma(zhat, n, nsigma, Pkk)

end

function sigmaPointsAlt(zhat::Vector{Float64}, n::Int, Pkk::Matrix{Float64})
    #zhat: state estimate
    #n: size of state vector
    #Wo: weighting of first sigma point
    #Pkk: covariance of state estimate
    #returns: chi (n x n*2+1) matrix with columns as sigma points

    chi = zeros(n, n * 2 + 1)

    Phalf = cholesky(Pkk);
    offset = sqrt(n);

    chi[:, 1] = zhat
    for i = 1:n
        
        chi[:, 1 + i] = zhat + offset * Phalf.L[:,i]
        chi[:, 1 + n + i] = zhat - offset * Phalf.L[:,i]

    end

    return chi

end

function sigmaPoints(zhat::Vector{Float64}, n::Int, Pkk::Matrix{Float64})

    return sigmaPoints(zhat, n, 1.0 - n/3, Pkk)

end

function sigmaPointsWeights(zhat::Vector{Float64}, n::Int, Pkk::Matrix{Float64})

    W0 = 1.0 - n/3
    return sigmaPoints(zhat, n, W0, Pkk), [W0, (1-W0)/(2*n)]

end

function ntheta2quatCov(n::Vector{Float64}, theta::Float64, dn::Vector{Float64}, dtheta::Float64)
    #n: rotation axis
    #theta: rotation angle
    #dn: diagonal elements of n covariance (assuming no cross covariance)
    #dtheta: covariance of angle theta
    #returns: dq diagonal of covariance of quaternion 

    #not sure how riggorus this is?

    dq_theta = vcat(0.5 * cos(theta/2) * n * dtheta, -0.5 * sin(theta/2) * dtheta)
    dq = dq_theta + vcat(dn, 0)

    return dq .* dq

end

#sensor function
function yhat(t::Float64, zhat::Vector{Float64}, simParam::sim)

    return yhat(t, zhat, stateDerivative(t, zhat, simParam))

end

function yhat(t::Float64, zhat::Vector{Float64}, dz::Vector{Float64})

    accel = dz[4:6]
    w_gyro = zhat[11:13]
    x3 = zhat[3]
    Be_B = rotateFrame(Be_I, zhat[7:10]) #direction of the earths magnetic field measured in body frame components

    return [accel; w_gyro; x3; Be_B]


end

function R(t::Float64, yhat::Vector{Float64}, lv::rocket)

    accel_cov_factor = 0.08
    accel_cov_base = 0.02
    gyro_cov_factor = 0.05
    baro_cov = 100
    mag_cov = 0.05 #measurement error of the magnetic field directions

    return diagm([abs(yhat[1]) * accel_cov_factor + accel_cov_base, abs(yhat[2]) * accel_cov_factor + accel_cov_base, abs(yhat[3]) * accel_cov_factor + accel_cov_base, abs(yhat[4]) * gyro_cov_factor, abs(yhat[5]) * gyro_cov_factor, abs(yhat[6]) * gyro_cov_factor, baro_cov, mag_cov,mag_cov, mag_cov])
    
end

#prediction step of unscented kalman filter
function ukf_step(tk::Float64, zkk::Vector{Float64}, Pkk::Matrix{Float64}, yk1::Vector{Float64}, Gw::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, dt::Float64, simParam::sim, nsigma::Float64)
    #tk: systetm time at known step (1x1)
    #zkk: system state at current time (known) (R_STATE_SIZEx1)
    #Pkk: system covariance (STATE_SIZE x STATE_SIZE)
    #yk1: measurement at time t + dt (Y_SIZE x 1)
    #Gw: matrix that transforms process noise covarince to state covariance (STATE_SIZE x W_SIZE)
    #Q: process noise covriance matrix (W_SIZE x W_SIZE)
    #R: Sensor noise covariance matrix (Y_SIZE x Y_SIZE)
    #dt: time step (1x1)

    updateMassState!(tk, simParam.rocket.massData, simParam.rocket.motorData) #update mass with current time
    dz(t,zi) = stateDerivative(t, zi, simParam) #find state derivative

    #calculating the weights to be used
    w0S, w0C, wi = nsigma_weights(nsigma)

    #PREDICTION

    qkk = zkk[7:10]

    ztildakk = z2ztilda(zkk, qkk)
    
    chitilda = sigmaPoints_nsigma(ztildakk, STATE_SIZE, nsigma, Pkk) #make sigma points based on covarince

    chi = zeros(R_STATE_SIZE, size(chitilda)[2])

    for i = 1:size(chi)[2]

        chi[:,i] = ztilda2z(chitilda[:,i], qkk)

    end

    #propgate each sigma point
    fchi = zeros(size(chi)[1], size(chi)[2])

    #sigma points for prediction step only
    zk1k = zeros(R_STATE_SIZE,1) #make state prediction to retrieve propgated q
    w = w0S
    for i = 1:size(chi)[2]

        fchi[:,i] = rk4Step(dz, tk, chi[:,i], dt) #propgate each sigma point

        if i == 2
            w = wi
        end

        #add new sigma point to state vector
        zk1k = zk1k + w * fchi[:,i]

    end
    #qk1k = zk1k[7:10] #find propogated quaternion
    qk1k = fchi[7:10,1]

    #transform each propogated sigma point back to estimator state
    chitildak1k = zeros(STATE_SIZE, size(chitilda)[2])

    ztildak1k = zeros(STATE_SIZE)
    w = w0S
    for i = 1:size(chi)[2]

        if i == 2
            w = wi
        end

        chitildak1k[:,i] = z2ztilda(fchi[:,i], qk1k)

        #add new sigma point to state vector
        ztildak1k = ztildak1k + w * chitildak1k[:,i]
        

    end

    #covariance predict
    Pk1k = Gw * Q * Gw' #initialize with process noise addition
    w = w0C #set weight to the first index
    for (index,fadd) in enumerate(eachcol(chitildak1k))

        #flip to other weight after the first index is added
        if index == 2
            w = wi
        end
        
        fadd = Hermitian((fadd - ztildak1k) * (fadd - ztildak1k)')
        Pk1k = Matrix(Hermitian(Pk1k + w * fadd))

    end

    #KALMAN UPDATE
    
    #redo sigma points with the predicted covariance
    chitildak1k = sigmaPoints_nsigma(ztildak1k, STATE_SIZE, nsigma, Pk1k) #new sigma points (already propgated as we are using ztildak1k)

    #transform new tilda sigma points to sigma points that can be plugged into yhat
    chik1k = zeros(R_STATE_SIZE, size(chitilda)[2])
    for i = 1:size(chik1k)[2]

        chik1k[:,i] = ztilda2z(chitildak1k[:,i], qk1k)

    end

    hchi = zeros(Y_SIZE, size(chik1k)[2])

    for i = 1:size(chik1k)[2]

        hchi[:,i] = yhat(tk+dt, chik1k[:,i], dz(tk+dt, chik1k[:,i])) #predicted sensor measurement for each propogated sigma point

    end

    #predicted sensor reading
    yk1k = zeros(Y_SIZE)
    w = w0S
    for (index,hadd) in enumerate(eachcol(hchi))

        #flip to other weight after the first index is added
        if index == 2
            w = wi
        end

        yk1k = yk1k + w * hadd

    end

    #yk1k = hchi[:,1]

    # print("yk1k: ")
    # println(yk1k)
    
    #covariance update
    Pk1k_yy = R
    Pk1k_zy = zeros(STATE_SIZE, Y_SIZE)
    w = w0C #set weight to the first index
    for i = 1:size(chi)[2]

        #flip to other weight after the first index is added
        if i == 2
            w = wi
        end

        Pk1k_yy = Pk1k_yy + w * ((hchi[:,i] - yk1k) * (hchi[:,i] - yk1k)')
        Pk1k_zy = Pk1k_zy + w * (chitildak1k[:,i] - ztildak1k) * (hchi[:,i] - yk1k)'
        # print(i)
        # print(": Pk1k_zy: ")
        # println(isposdef(Pk1k_zy))

    end
 
    # print("yk1-yk1k")
    # println(yk1 - yk1k)

    ztildak1k1 = ztildak1k + Pk1k_zy * (Pk1k_yy \ (yk1 - yk1k))
    zk1k1 = ztilda2z(ztildak1k1, qk1k)
    Pkadjust = Hermitian(Pk1k_zy * (Pk1k_yy \ (Pk1k_zy')))

    factor = 1.0
    # print("Pkadjust: ")
    # println(isposdef(Pkadjust))
    Pk1k1 = Pk1k - Pkadjust

    while(!(isposdef(Pk1k1)) && factor > 0.01)
        factor = factor * 0.9999
        print("Factor: ")
        println(factor)
        Pk1k1 = Pk1k - factor * Pkadjust
    end


    return zk1k1, Pk1k1

end

function ztilda2z(ztilda::Vector{Float64}, qbase::Vector{Float64})
    #ztilda: estimator state 
    #qbase: quaternion that ds is based from 
    #returns: z (system state) with updated quaternion 

    a = 0.5
    f = 2 * (a + 1)

    dq = rev2quatError(ztilda[7:9], f, a)
    qnew = quatProd(dq, qbase)
    

    return [ztilda[1:6]; qnew; ztilda[10:12]]

end
function z2ztilda(z::Vector{Float64}, qbase::Vector{Float64})
    #z: System state vector (R_STATE_SIZE x 1)
    #qbase: 
    #returns: ztilda --> estimation state vector with ds = 0 (STATE_SIZE x 1)

    a = 0.5
    f = 2 * (a + 1)

    dq = quatProd(z[7:10], quatInv(qbase))
    ds = quatError2rev(dq, f, a)

    return [z[1:6]; ds; z[11:13]]

end

function ukf(simParam::sim, y::Matrix{Float64}, P0::Matrix{Float64}, nsigma::Float64)
    #simParam: sim struct that describes the expected behavior of the system (expected wind + ideal motor)
    #y: matrix containing all the generated sensor data (length(tspan) x Y_SIZE)
    #P0: initial covariance matrix
    #nsigma: parameter for chosing sigma points --> distance of sigma points from mean

    num_steps = length(simParam.simInputs.tspan)
    tspan = simParam.simInputs.tspan
    dt = tspan[2] - tspan[1] #assumes equal time spacing


    zhat = zeros(R_STATE_SIZE, num_steps) #system state
    Phat = zeros(STATE_SIZE, STATE_SIZE, num_steps)
    
    zhat[:,1] = simParam.simInputs.z0
    Phat[:,:,1] = P0

    for k = 1:num_steps - 1

        print("STEP: ")
        println(k)
        

        Gw = getGw(tspan[k], zhat[:,k], dt, simParam)
        Q = getQ(tspan[k], zhat[:,k], simParam.rocket)
        Rw = R(tspan[k], yhat(tspan[k], zhat[:,k], simParam), simParam.rocket)
        zhat[:,k+1], Phat[:,:,k+1] = ukf_step(tspan[k], zhat[:,k], Phat[:,:,k], y[k+1,:], Gw, Q, Rw, dt, simParam, nsigma)

    end

    return zhat, Phat

end

function testDataRun(simParam::sim, wind::windData, thrustVar::Float64)

    simParamProcessNoise = copy(simParam)
    simParamProcessNoise.simInputs.windData = wind
    simParamProcessNoise.simInputs.thrustVar = thrustVar

    tspan, z = run(simParamProcessNoise)
    data = noisySensorData(tspan, z, simParamProcessNoise)

    return tspan, z, data
end

function noisySensorData(tspan::Vector{Float64}, z::Matrix{Float64}, simParam::sim)

    data = zeros(size(z)[1], Y_SIZE)
    for i = 1:size(data)[1]

        data[i, :] = yhat(tspan[i], z[i,:], simParam)'
        Rk = R(tspan[i], data[i,:], simParam.rocket)
        
        if(isposdef(Rk))
            dist = MvNormal(zeros(Y_SIZE), Rk)
            noise = rand(dist)
            data[i, :] = data[i, :] + noise
        end

    end

    return data

end

function printMat(a::Matrix{Float64})

    m,n = size(a)

    for i = 1:m
        for j = 1:n
            @printf("%.4f", a[i,j])
            print(", ")
        end
        println()
    end

end

function plotEstimator(tspan::Vector{Float64}, ztrue::Vector{Float64}, zhat::Vector{Float64}, Pu::Vector{Float64}, t::String)

    pygui(true)
    upperBound = zhat + 2 * sqrt.(Pu)
    lowerBound = zhat - 2 * sqrt.(Pu)

    figure()
    title(t)
    plot(tspan, ztrue)
    plot(tspan, zhat)
    plot(tspan, upperBound, "--")
    plot(tspan, lowerBound, "--")
    legend(["ztrue", "zhat"])

end


function plotEstimatorError(tspan::Vector{Float64}, ztrue::Vector{Float64}, zhat::Vector{Float64}, Pu::Vector{Float64}, t::String)

    pygui(true)
    upperBound = 2 * sqrt.(Pu)
    lowerBound = -2 * sqrt.(Pu)

    figure()
    title(t)
    plot(tspan, (ztrue - zhat))
    axhline(y = 0.0, color="red")
    plot(tspan, upperBound, "--")
    plot(tspan, lowerBound, "--")
    legend(["error", "ztrue"])

end

let 
    winds = [0.0 0.0 0.0 0; 
             0 0  0 0;
             0  0  0 0]
    h = [0.0, 1000, 2000, 3000]

    trueWind = windData(h, winds)
    nsigma = 1.0

    #expected data
    simRead = readJSONParam("simParam.JSON")
    setWindData!(simRead.simInputs, [0.0, 1000.0],  [0.0 0.0; 0.0 0.0; 0.0 0.0])

    z0 = simRead.simInputs.z0

    ds = [0.02, 0.02, 0.02]
    dx = [10.0,10.0,10.0]
    dv = [0.01,0.01,0.01]
    dw = [0.01,0.01,0.01]

    P0 = diagm(vcat(dx,dv,ds,dw))

    


    # tspan, expected_z = run(simRead) 

    tspan, ztrue, y = testDataRun(simRead, trueWind, 1.0)

    # getQuiverPlot_py(expected_z, 1)

    getQuiverPlot_py(ztrue, 1)

    zhat, Phat = @timev ukf(simRead, y, P0, nsigma)

    getQuiverPlot_py(transpose(zhat), 1)

    j = 3
    plotEstimatorError(tspan, ztrue[:,j], zhat[j,:], Phat[j,j,:], "x3")

    j = 6
    plotEstimatorError(tspan, ztrue[:,j], zhat[j,:], Phat[j,j,:], "v3")

    j = 4
    plotEstimatorError(tspan, ztrue[:,j], zhat[j,:], Phat[j,j,:], "v1")
    




    ##old tests

    # q0 = z0[7:10]

    # ztildakk = z2ztilda(z0, q0)
    
    # chitilda = sigmaPoints_nsigma(ztildakk, STATE_SIZE, nsigma, P0) #make sigma points based on covarince

    # chi = zeros(R_STATE_SIZE, size(chitilda)[2])

    # for i = 1:size(chi)[2]

    #     chi[:,i] = ztilda2z(chitilda[:,i], q0)

    # end

    # b3Imat = zeros(3, size(chi)[2])

    # for i = 1:size(chi)[2]

    #     b3Imat[:,i] = rotateFrame([0,0,1.0], quatInv(chi[7:10,i]))

    # end

    # pygui(true)
    # scatter3D(b3Imat[1,:],b3Imat[2,:],b3Imat[3,:])
    # scatter3D(0,0,0)

        #needs to match that in the sim param file
    # n0 = [0,1.0,0]
    # theta0 = 5 * pi/180

    #code for displaying whole matrix thing 
    #show(IOContext(stdout, :limit=>false), MIME"text/plain"(), <MATRIX>)

    # z_test = [0,0,1500,10,10,100.0,0,0,0,1,0,0,0]
    # Pkk = diagm([5, 5, 5, 5, 5, 5, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.5])
    # dt = 0.05

    #testing sigmaPoints generator 

    #chi = sigmaPoints(z_test, STATE_SIZE, Pkk)


    #testing the Gw function
    # Gw_test = getGw(0.5, z_test, .05, rocket)
    # println(Gw_test[4:6,:])

end