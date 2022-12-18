#Sigma Point Filter (UKF) Implimentation for Cornell Rocketry Real-Time State Estimation
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
#ds: generalized Rodrigues error vector representing error from the quaternion in the system state
#w:IwB rotation of B frame wrt to the inertial frame 11,12,13


###Estimator#####

#y: Measurement 11x1: [y_accel, y_gyro, y_barox3, y_magb]
#R: Sensor noise covariance 10x10
#w: [ww1, ww2, ww3, wt, wF1, wF2, wF3, wm1, wm2, wm3] 
#Q: Process noise covariance WSIZExWSIZE

include("6dof.jl")

const global STATE_SIZE = 13 #estimator state size
const global RE_STATE_SIZE = 14 #estimator and rocket state (rocket state + parameters)
const global W_SIZE = 11;
const global Y_SIZE = 10;
const global Be_I = [0.0;1.0;0.0] #direction of the Earths magnetic field in the inertial frame
global a = 1.0 #constant used in the back and forth between generalizedd Rodrigues error vector and error quaternion


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

    Gw[13, 11] = 1.0;  #process noise for the motor thrust parameter
    
    return Gw

end

function getQ(t::Float64, z::Vector{Float64}, lv::rocket)

    Qwind = getWindVar(t, z[3])
    thrust = motorThrustMass(t, lv.motorData, lv.massData.initalMass[2])[1]
    Ig_B = getIg(lv.massData)
    thrustVariation = sqrt(.16)

    addativeForceProcessNoise = 500.0 #Newtons^2
    addativeMomentProcessNoise = 100.0 #(Nm)^2
    thrustParamCov = 0.04

    Q = zeros(W_SIZE,W_SIZE)
    Q[1:3,1:3] = Qwind
    Q[4,4] = (thrustVariation * thrust)^2
    Q[5:7,5:7] = addativeForceProcessNoise * diagm(ones(3))
    Q[8:10,8:10] = addativeMomentProcessNoise * diagm(ones(3))
    Q[11,11] = thrustParamCov

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

    accel_cov_factor = 0.0005 #(m^2s^-4)
    gyro_cov_factor = 0.0005 #(s^-2)
    baro_cov = 9 #(m^2)
    mag_cov = 0.01 #measurement error of the magnetic field directions

    rvec = [ones(3) * accel_cov_factor; ones(3) * gyro_cov_factor; baro_cov; ones(3) * mag_cov]

    return diagm(rvec)
    
end

#prediction step of unscented kalman filter
function ukf_step(tk::Float64, zkk::Vector{Float64}, Pkk::Matrix{Float64}, yk1::Vector{Float64}, Gw::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, dt::Float64, simParam::sim, nsigma::Float64)
    #tk: systetm time at known step (1x1)
    #zkk: system state at current time (known) (RE_STATE_SIZEx1)
    #Pkk: system covariance (STATE_SIZE x STATE_SIZE)
    #yk1: measurement at time t + dt (Y_SIZE x 1)
    #Gw: matrix that transforms process noise covarince to state covariance (STATE_SIZE x W_SIZE)
    #Q: process noise covriance matrix (W_SIZE x W_SIZE)
    #R: Sensor noise covariance matrix (Y_SIZE x Y_SIZE)
    #dt: time step (1x1)

    dz(t,zi) = stateDerivative(t, zi, simParam) #find state derivative

    #calculating the weights to be used
    w0S, w0C, wi = nsigma_weights(nsigma)

    #PREDICTION

    qkk = zkk[7:10]
    # simParam.simInputs.thrustVar = zkk[14] #set thrustVar as alpha
    # print("zkk14: ")
    # println(zkk[14])

    # print("ThrustVarStartUKF: ")
    # println(simParam.simInputs.thrustVar)

    ztildakk = z2ztilda(zkk, qkk)
    
    chitilda = sigmaPoints_nsigma(ztildakk, STATE_SIZE, nsigma, Pkk) #make sigma points based on covarince

    chi = zeros(RE_STATE_SIZE, size(chitilda)[2])

    for i = 1:size(chi)[2]

        chi[:,i] = ztilda2z(chitilda[:,i], qkk)

    end

    #propgate each sigma point
    fchi = zeros(size(chi)[1], size(chi)[2])

    #sigma points for prediction step only
    w = w0S
    for i = 1:size(chi)[2]

        simParam.simInputs.thrustVar = chi[14,i] #set thrustVar to parameter in each sigma-point
        fchi[1:R_STATE_SIZE,i] = rk4Step(dz, tk, chi[1:R_STATE_SIZE,i], dt) #propgate each sigma point
        fchi[R_STATE_SIZE+1:RE_STATE_SIZE, i] = chi[R_STATE_SIZE+1:RE_STATE_SIZE,i] #propgating the parameters is just copying them

    end
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
    chik1k = zeros(RE_STATE_SIZE, size(chitilda)[2])
    
    #chik1k = fchi
    for i = 1:size(chik1k)[2]

        chik1k[:,i] = ztilda2z(chitildak1k[:,i], qk1k)

    end

    updateMassState!(tk + dt, simParam.rocket.massData, simParam.rocket.motorData) #update mass with current time (tk + dt)

    hchi = zeros(Y_SIZE, size(chik1k)[2])

    for i = 1:size(chik1k)[2]

        simParam.simInputs.thrustVar = chik1k[14,i] #set thrustVar to sigma point value before calculating sensor value
        hchi[:,i] = yhat(tk+dt, chik1k[1:R_STATE_SIZE,i], dz(tk+dt, chik1k[1:R_STATE_SIZE,i])) #predicted sensor measurement for each propogated sigma point

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

    #covariance update
    Pk1k_yy = R
    Pk1k_zy = zeros(STATE_SIZE, Y_SIZE)
    w = w0C #set weight to the first index
    for i = 1:size(hchi)[2]

        #flip to other weight after the first index is added
        if i == 2
            w = wi
        end

        Pk1k_yy = Pk1k_yy + w * ((hchi[:,i] - yk1k) * (hchi[:,i] - yk1k)')
        Pk1k_zy = Pk1k_zy + w * (chitildak1k[:,i] - ztildak1k) * (hchi[:,i] - yk1k)'

    end

    ztildak1k1 = ztildak1k + Pk1k_zy * (Pk1k_yy \ (yk1 - yk1k))
    zk1k1 = ztilda2z(ztildak1k1, qk1k)
    Pkadjust = Hermitian(Pk1k_zy * (Pk1k_yy \ (Pk1k_zy')))
    Pk1k1 = Pk1k - Pkadjust

    # factor = 1.0
    # while(!(isposdef(Pk1k1)) && factor > 0.01)
    #     factor = factor * 0.9999999
    #     print("Factor: ")
    #     println(factor)
    #     Pk1k1 = Pk1k - factor * Pkadjust
    # end

    simParam.simInputs.thrustVar = zk1k1[14]

    return zk1k1, Pk1k1

end

function ztilda2z(ztilda::Vector{Float64}, qbase::Vector{Float64})
    #ztilda: estimator state (STATE_SIZE x 1)
    #qbase: quaternion that ds is based from 
    #returns: z (system state) with updated quaternion (RE_STATE_SIZE)

    f = 2 * (a + 1)

    dq = rev2quatError(ztilda[7:9], f, a)
    qnew = quatProd(dq, qbase)
    

    return [ztilda[1:6]; qnew; ztilda[10:STATE_SIZE]]

end

function z2ztilda(z::Vector{Float64}, qbase::Vector{Float64})
    #z: System state vector (RE_STATE_SIZE x 1)
    #qbase: 
    #returns: ztilda --> estimation state vector with ds = 0 (STATE_SIZE x 1)

    f = 2 * (a + 1)

    dq = quatProd(z[7:10], quatInv(qbase))
    ds = quatError2rev(dq, f, a)

    return [z[1:6]; ds; z[11:RE_STATE_SIZE]]

end

function ukf(simParam::sim, y::Matrix{Float64}, z0::Vector{Float64}, P0::Matrix{Float64}, nsigma::Float64)
    #simParam: sim struct that describes the expected behavior of the system (expected wind + ideal motor)
    #y: matrix containing all the generated sensor data (length(tspan) x Y_SIZE)
    #z0: starting value for state estimator
    #P0: initial covariance matrix
    #nsigma: parameter for chosing sigma points --> distance of sigma points from mean

    num_steps = length(simParam.simInputs.tspan)
    tspan = simParam.simInputs.tspan
    dt = tspan[2] - tspan[1] #assumes equal time spacing


    zhat = zeros(RE_STATE_SIZE, num_steps) #system state
    Phat = zeros(STATE_SIZE, STATE_SIZE, num_steps)
    
    zhat[:,1] = z0
    Phat[:,:,1] = P0

    for k = 1:num_steps - 1

        print("STEP: ")
        println(k)
        print("Thrust Var")
        println(zhat[14,k])
        
        updateMassState!(tspan[k], simParam.rocket.massData, simParam.rocket.motorData)


        Gw = getGw(tspan[k], zhat[:,k], dt, simParam)
        Q = getQ(tspan[k], zhat[:,k], simParam.rocket)
        Rw = R(tspan[k+1], yhat(tspan[k+1], zhat[:,k+1], simParam), simParam.rocket)
        zhat[:,k+1], Phat[:,:,k+1] = ukf_step(tspan[k], zhat[:,k], Phat[:,:,k], y[k+1,:], Gw, Q, Rw, dt, simParam, nsigma)
        

    end

    return zhat, Phat

end

function testDataRun!(simParam::sim, wind::windData, thrustVar::Float64)

    simParam.simInputs.windData = wind
    simParam.simInputs.thrustVar = thrustVar
    
    return testDataRun!(simParam)
    
end

function testDataRun!(simParam::sim)

    tspan, z = run(simParam)
    data = noisySensorData(tspan, z, simParam)

    return tspan, z, data

end

function noisySensorData(tspan::Vector{Float64}, z::Matrix{Float64}, simParam::sim)

    data = zeros(size(z)[1], Y_SIZE)
    for i = 1:size(data)[1]

        updateMassState!(tspan[i], simParam.rocket.massData, simParam.rocket.motorData)
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
    plot(tspan, upperBound, color="green", "--")
    plot(tspan, lowerBound, color="green", "--")
    legend(["error", "ztrue", "2σ bound"])
    xlabel("time (s)")
    ylabel(t)

end

function plotEstimator(tspan::Vector{Float64}, ztrue::Vector{Float64}, zhat::Vector{Float64}, t::String)

    pygui(true)
    figure()
    title(t)
    plot(tspan, ztrue)
    plot(tspan, zhat)
    legend(["ztrue", "zhat"])
    xlabel("time (s)")
    ylabel(t)

end


function plotEstimatorError(tspan::Vector{Float64}, ztrue::Vector{Float64}, zhat::Vector{Float64}, Pu::Vector{Float64}, t::String)

    pygui(true)
    upperBound = 2 * sqrt.(Pu)
    lowerBound = -2 * sqrt.(Pu)

    figure()
    title(t)
    plot(tspan, (ztrue - zhat))
    axhline(y = 0.0, color="red")
    axvline(x = 2.483, color="gray")
    plot(tspan, upperBound, color="green", "--")
    plot(tspan, lowerBound, color="green",  "--")
    legend(["error", "ztrue", "burnout", "2σ bound"])
    xlabel("time (s)")
    ylabel(t * "Error")

end


let 
    truewinds = [0.0 5.0 0.0 0; 
             0 0  0 0;
             0  0  0 0]
    expectedwinds = [0.0 5.0 0.0 0; 
                0 0  0 0;
                0  0  0 0]

    h = [0.0, 1000, 2000, 3000]

    trueWind = windData(h, truewinds)
    expectedWind = windData(h, expectedwinds)
    

    #expected data
    simTrue = readJSONParam("simParam.JSON") #simParm used for generating data and true values
    simExpected = readJSONParam("simParam.JSON") #simParam used for the estimator 

    simExpected.simInputs.windData = expectedWind

    #setting parameters for data generation
    thrustVarMean = 1.05 #mean multiple of ideal thrust at each time step
    thrustVarCov = .0025 #variation of multiple of ideal thrust at each time step
    simTrue.simInputs.windData = trueWind
    simTrue.simInputs.isDataGen = true
    simTrue.simInputs.dataGenThrustVar = rand(Normal(thrustVarMean, thrustVarCov), length(simTrue.simInputs.tspan))

    #generate ztrue and data with true wind and thrust variation
    tspan, ztrue, y = testDataRun!(simTrue)

    getQuiverPlot_py(ztrue, 1)

    #setting up initial covariance
    ds = [0.02, 0.02, 0.02]
    dx = [10.0,10.0,5.0]
    dv = [0.01,0.01,0.01]
    dw = [0.01,0.01,0.01]
    dTv = 1e-4
    P0 = diagm(vcat(dx,dv,ds,dw, dTv))

    z0 = simExpected.simInputs.z0
    initalThrustVarEstimate = 1.0 #you assume it is going to perform nominally
    simExpected.simInputs.thrustVar = initalThrustVarEstimate  #set thrust var to be correct inital value

    #setting up z0 vector with the parameters at the end
    z0 = vcat(z0, initalThrustVarEstimate) 
    println(simExpected.simInputs.isDataGen)

    nsigma = 1.0
    zhat, Phat = @timev ukf(simExpected, y, z0, P0, nsigma)

    getQuiverPlot_py(transpose(zhat), 1)

    j = 3
    plotEstimatorError(tspan, ztrue[:,j], zhat[j,:], Phat[j,j,:], "x3")

    j = 6
    plotEstimatorError(tspan, ztrue[:,j], zhat[j,:], Phat[j,j,:], "v3")

    j = 11
    plotEstimatorError(tspan, ztrue[:,j], zhat[j,:], Phat[j,j,:], "w1")

    plotnumber = 235
    plotEstimatorError(tspan[1:plotnumber], ones(plotnumber) * thrustVarMean, zhat[14,1:plotnumber], Phat[13,13,1:plotnumber], "thrustVar")
    




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