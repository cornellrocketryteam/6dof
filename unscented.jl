#Unscented Kalman Filter Implimentation for Cornell Rocketry Real-Time State Estimation
#Developed by Sam Noles in 2022

#PyPlot Plotting Examples: https://gist.github.com/gizmaa/7214002#multiaxis

#########info##############

#system state: z = [r v q w]
#r: position of R (tip of nose cone) wrt to O (ground) in inertial frame 1,2,3
#v: inertial frame derivative of r in intertial frame coordinates 4,5,6
#q: quaternion representing rotation from inertial to body frame 7,8,9,10
    #aka apply q to a vector written in I cordinates to get the same vector written in body frame coordinates
#w:IwB rotatoin of B frame wrt to the inertial frame 11,12,13

###Estimator#####

#y: Measurement 11x1: [y_accel, y_gyro, y_barox3, y_magq]
#R: Sensor noise covariance 11x11
#w: [ww1, ww2, ww3, wt]
#Q: Process noise covariance 4x4

include("6dof.jl")

const global W_SIZE = 4;
const global Y_SIZE = 11;


#functions 

#implimentation of Gw(t,z)
function getGw(t::Float64, z::Vector{Float64}, dt::Float64, rocket::rocket)
    #t: current time
    #z: state vector
    #dt: time step
    #rocket: rocket object that has all the rocket info
    #return: Gw STATE_SIZE x 4 matrix w = [ww1, ww2, ww3, wt]

    Gw = zeros(STATE_SIZE, W_SIZE);
    m = sum(rocket.massData.currentMass)

    #contribution from thrust uncertainty
    Gw[6, 4] = dt / m

    vAOI_I = getWind(t, z[3])[1]
    vRAI_I = getVRA(z[4:6], vAOI_I)
    mag_vRAI_I = norm(vRAI_I)
    aoa = calcAoA(vRAI_I, z[7:10])
    mach = calcMach(vRAI_I, z[3])


    #wind uncertainty
    dfd_dwd = dt * 0.5 * getCd(aoa, mach, rocket.aeroData) * getA(aoa, rocket.aeroData) * expAtm(z[3]) * (vRAI_I * transpose(vRAI_I) / mag_vRAI_I + (I * mag_vRAI_I))/m

    Gw[4:6, 1:3] = dfd_dwd;
    
    return Gw

end

function getQ(t::Float64, z::Vector{Float64}, lv::rocket)

    Qwind = getWind(t, z[3])[2]
    thrust = motorThrustMass(t, lv.motorData, lv.massData.initalMass[2])[1]
    thrustVariation = getThrustVar(t)

    Q = zeros(4,4)
    Q[1:3,1:3] = Qwind
    Q[4,4] = thrustVariation * thrust

    return Q
end


function sigmaPoints(zhat::Vector{Float64}, n::Int, Wo::Float64, Pkk::Matrix{Float64})
    #zhat: state estimate
    #n: size of state vector
    #Wo: weighting of first sigma point
    #Pkk: covariance of state estimate
    #returns: chi (n x n*2+1) matrix with columns as sigma points

    chi = zeros(n, n * 2 + 1)
    Phalf = cholesky(Pkk);
    offset = sqrt(n/(1-Wo))

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
function yhat(t::Float64, zhat::Vector{Float64}, lv::rocket)

    return yhat(t, zhat, stateDerivative(t, zhat, lv))

end

function yhat(t::Float64, zhat::Vector{Float64}, dz::Vector{Float64})

    accel = dz[4:6]
    w_gyro = zhat[11:13]
    x3 = zhat[3]
    q = zhat[7:10]

    return [accel; w_gyro; x3; q]


end

function R(t::Float64, yhat::Vector{Float64}, lv::rocket)

    accel_cov_factor = 0.03
    gyro_cov_factor = 0.03
    baro_cov = 5
    mag_cov = 0.01 # no idea how to do this covariance 

    return diagm([yhat[1] * accel_cov_factor, yhat[2] * accel_cov_factor, yhat[3] * accel_cov_factor, yhat[4] * gyro_cov_factor, yhat[5] * gyro_cov_factor, yhat[6] * gyro_cov_factor, baro_cov, mag_cov,mag_cov, mag_cov, mag_cov])
    
end

#prediction step of unscented kalman filter
function ukf_step(tk::Float64, zkk::Vector{Float64}, Pkk::Matrix{Float64}, yk1::Vector{Float64}, Gw::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, dt::Float64, lv::rocket)
    #tk: systetm time at known step 
    #zkk: system state at current time (known)
    #Pkk: system covariance
    #yk1: measurement at time t + dt
    #Gw: matrix that transforms process noise covarince to state covariance
    #Q: process noise covriance matrix
    #R: Sensor noise covariance matrix
    #dt: time step

    dz(t,zi) = stateDerivative(t, zi, lv.aeroData, lv.massData, lv.motorData, lv.latLong)
    
    chi, W = sigmaPointsWeights(zkk, length(zkk), Pkk)

    #calculating states for uncented sums 
    fchi = zeros(size(chi)[1], size(chi)[2])
    hchi = zeros(Y_SIZE, size(chi)[2])



    for i = 1:size(chi)[2]

        fchi[:,i] = rk4Step(dz, tk, chi[:,i], dt)
        hchi[:,i] = yhat(tk, chi[:,i], dz(tk, chi[:,i]))

    end

    #PREDICTION
    zk1k = fchi[:,1]
    Pk1k = Gw * Q * Gw'
    wi = W[1] #set weight to the first index
    for (index,fadd) in enumerate(eachcol(fchi))

        #flip to other weight after the first index is added
        if index == 2
            wi = W[2]
        end

        Pk1k = Pk1k + wi * (fadd - zk1k) * transpose(fadd - zk1k)

    end

    #KALMAN UPDATE

    yk1k = hchi[:,1] #predicted sensor reading
    Sk1 = R
    C = zeros(STATE_SIZE, Y_SIZE)

    wi = W[1] #set weight to the first index
    for i = 1:size(chi)[2]

        #flip to other weight after the first index is added
        if index == 2
            wi = W[2]
        end

        Sk1 = Sk1 + wi * (hchi[:,i] - yk1k) * transpose(hchi[:,i] - yk1k)
        C = C + wi * (fchi[:,i] - zk1k) * transpose(hchi[:,i] - yk1k)

    end

    Sk1_factorized = cholesky(Sk1)
    zk1k1 = zk1k + C * Sk1_factorized \ (yk1 - yk1k)
    Pk1k1 = Pk1k - C * Sk1_factorized \ transpose(C)

    return zk1k1, Pk1k1

end


let 

    simRead = readJSONParam("simParam.JSON")
    rocket = simRead.rocket

    z_test = [0,0,1500,10,10,100.0,0,0,0,1,0,0,0]
    Pkk = diagm([5, 5, 5, 5, 5, 5, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.5])
    dt = 0.05

    #testing kalman step

    ti = 1.0

    ukf_step(ti, z_test, Pkk, getGw(ti, z_test, dt, rocket), getQ(ti, z_test, rocket), R(ti, h(ti,z_test), rocket), dt, rocket)


    #testing sigmaPoints generator 

    #chi = sigmaPoints(z_test, STATE_SIZE, Pkk)


    #testing the Gw function
    # Gw_test = getGw(0.5, z_test, .05, rocket)
    # println(Gw_test[4:6,:])

end