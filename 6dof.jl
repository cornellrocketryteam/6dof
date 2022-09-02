#6DoF modoel for the Cornell Rocketry Team 
#Written in 2022 by Samuel Noles

#########info##############

#system state: z = [r v q w]
#r: position of R (tip of nose cone) wrt to O (ground) in inertial frame 1,2,3
#v: inertial frame derivative of r in intertial frame coordinates 4,5,6
#q: quaternion representing rotation from inertial to body frame 7,8,9,10
    #aka apply q to a vector written in I cordinates to get the same vector written in body frame coordinates
#w:IwB rotatoin of B frame wrt to the inertial frame 11,12,13

#system time 
#t: t = 0 at motor ignition

using LinearAlgebra
using Interpolations
using StructArrays
using PyPlot
using Optim


########structs################

struct aeroDataPoint
    #independent variables
    AoA::Float64 #angle of attack
    Mach::Float64 #Mach numbers

    #dependent variables
    Cd::Float64 #coefficient of drag
    Cl::Float64 #coefficient of lift
    COP::Float64 #center of pressure
    A::Float64 #cross sectional area
end

mutable struct massElement
    
    comPos::Vector{Float64} #postition of center of mass in body frame coordinates (should be negative)
    initalMass::Float64 #total inital mass of component
    currentMass::Float64 #updated mass of component (propellant as rocket burns)
    I::Matrix{Float64} #3x3 body frame inertia tensor
    Icalc::Function #function used to update I

end

struct windPoint
    height::Float64
    windVector::Array{Float64}
end

mutable struct aeroCharacterization
    AoA::Vector{Float64}
    Mach::Vector{Float64}

    Cd::Matrix{Float64}
    Cl::Matrix{Float64}
    COP::Matrix{Float64}
    A::Vector{Float64}

    aeroCharacterization() = new([],[],Array{Float64}(undef, 0, 0),Array{Float64}(undef, 0, 0),Array{Float64}(undef, 0, 0),[])

end

#########functions (majority of project)#########

#mutating function that adds data to aeroCharacterization data set
function addData!(previous::aeroCharacterization, dataPoint::aeroDataPoint)

    #previous: !modifies! aeroDataCharacterization to describe a body
    #dataPoint: aeroDatapoint to add to the aeroDataCharacterization

    tol = 1e-10
    aoaIndex = getIndexWithTol(previous.AoA, dataPoint.AoA, tol)
    machIndex = getIndexWithTol(previous.Mach, dataPoint.Mach, tol)

    Cd = previous.Cd

    if(aoaIndex < 0)
        aoaIndex = aoaIndex * -1
        previous.AoA = [previous.AoA[1:aoaIndex-1]; dataPoint.AoA; previous.AoA[aoaIndex:end]]
        previous.Cd = addRow(previous.Cd, aoaIndex)
        previous.Cl = addRow(previous.Cl, aoaIndex)
        previous.COP = addRow(previous.COP, aoaIndex)
        previous.A = [previous.A[1:aoaIndex-1]; 0; previous.A[aoaIndex:end]]
    end

    if(machIndex < 0)
        machIndex = machIndex * -1
        previous.Mach = [previous.Mach[1:machIndex-1]; dataPoint.Mach; previous.Mach[machIndex:end]]
        previous.Cd = addCol(previous.Cd, machIndex)
        previous.Cl = addCol(previous.Cl, machIndex)
        previous.COP = addCol(previous.COP, machIndex)
    end

    previous.Cd[aoaIndex,machIndex] = dataPoint.Cd
    previous.Cl[aoaIndex,machIndex] = dataPoint.Cl
    previous.COP[aoaIndex,machIndex] = dataPoint.COP
    previous.A[aoaIndex] = dataPoint.A

end

#interpolate A from aeroCharacterization
function getA(AoA::Float64, data::aeroCharacterization)
    #AoA: input angle of attack
    #data: struct array of aeroDataPoints
    #return: interpolated area

    interp = LinearInterpolation(data.AoA, data.A)
    return interp(AoA);

end

#interpolate Cd from aeroCharacterization
function getCd(AoA::Float64, Mach::Float64, data::aeroCharacterization)
    #AoA: input angle of attack
    #Mach: input mach numbersString
    #data: struct array of aeroDataPoints
    #return: interpolated Cd

    interp = LinearInterpolation((data.AoA, data.Mach), data.Cd)
    return interp(AoA, Mach)

end

#interpolate Cl from aeroCharacterization
function getCl(AoA::Float64, Mach::Float64, data::aeroCharacterization)
    #AoA: input angle of attack
    #Mach: input mach numbersString
    #data: struct array of aeroDataPoints
    #return: interpolated Cl

    interp = LinearInterpolation((data.AoA, data.Mach), data.Cl)
    return interp(AoA, Mach)

end

#interpolate COP from aeroCharacterization
function getCOP(AoA::Float64, Mach::Float64, data::aeroCharacterization)
    #AoA: input angle of attack
    #Mach: input mach numbersString
    #data: struct array of aeroDataPoints
    #return: interpolated COP

    interp = LinearInterpolation((data.AoA, data.Mach), data.COP)
    return [0;0;interp(AoA, Mach)]

end
    
#finds index of value in array or index of where value would be if inserted
function getIndexWithTol(array::Vector, value::Float64, tol::Float64)
    #array: array (assumed to be sorted)
    #value: value being searched for
    #returns: index of first array element that is within tol of the given value
        # negative index implies value does not already exist but the value would exist at that index if added

    for (index, val) in enumerate(array)
        if abs(value - val) < tol
            return index
        end
    end

    #not already in array, find where it will go given array is sorted
    lowIndex = 1
    while lowIndex <= length(array) && value > array[lowIndex]
        lowIndex = lowIndex + 1
    end
    return -lowIndex
    
end

#reads in motor data and returns array of time in first column and thrust in the second
function readMotorData(thrustCurveFilePath::String)
    # thrustCurveFileName: file path to .eng file with thrust curve data
    ## reads in .eng files, modifies motorData to fill it with time/thrust data
    motorData = [0,0]'

    file = open(thrustCurveFilePath, "r")

    headerLine = readline(file)

    while !eof(file)
        line = readline(file)

        ##remove blank space at beginning
        while line[1] == ' '
            line = SubString(line, 2)
        end

        numbersString = split(line, ' ')
        numbers = [parse(Float64, numbersString[1]), parse(Float64, numbersString[2])]'
        motorData = vcat(motorData, numbers)
    end
    return motorData
end

function readRRC3Data(flightDataFilePath::String)

    mpf = 0.3048 #meters per foot conversion

    flightData = [0,0]'

    file = open(flightDataFilePath, "r")

    headerline = readline(file)

    while !eof(file)
        line = readline(file)
        
        numbersString = split(line, ',')
        dataLine = [parse(Float64, numbersString[1]), parse(Float64, numbersString[2]) * mpf]'
        flightData = vcat(flightData, dataLine)
    end

    flightData = flightData[2:end, :]
    return flightData
end

#gives mass and thrust of motor for given time after ignition
function motorThrustMass(t::Float64, motorData::Matrix{Float64}, initalPropMass::Float64)

    #t: time since motor ignition (s)
    #motorData: matrix of motor thrust time data
    #initalPropMass: mass of propellant in fully loaded motor
    #return: thrust (N), mass of propellant remaining
    ## unis must match, between data file and other


    #get motor burn time

    burnTime = motorData[size(motorData)[1], 1]

    if t >= burnTime
        return 0,0
    end

    #find range that t is in of the given datapoints in motorData

    indexLow = 1
    while (indexLow < (size(motorData)[1] - 1) && t > motorData[indexLow+1,1])
        indexLow = indexLow + 1
    end

    return motorData[indexLow,2] + (motorData[indexLow + 1, 2] - motorData[indexLow,2])/(motorData[indexLow + 1, 1] - motorData[indexLow,1]) * (t - motorData[indexLow, 1]), initalPropMass - initalPropMass * t / burnTime
end

#atmosphereic density as a function of height above sea level
function expAtm(h::Float64)

    #h: height above sea level (m)
    #returns ρ: density of atmosphere at that height (kg/m^3)

    ρ(ρ0, h, h0, H) = ρ0 * exp((h0-(h/1000))/H)

    return ρ(1.225, h, 0, 7.249) #Currently written for values of h between 0-25km add lookup table for more

end

##atmosphereic density as a function of position 
function expAtm_r(rRO_I)

    rOC_I = [0.0;0.0; 6.3781e6]
    rRC_I = rOC_I + rRO_I;
    h = norm(rRC_I) - norm(rOC_I)

    return expAtm(h)

end

#integrator step function
function rk4Step(dz::Function, ti::Float64, zi::Vector{Float64}, dt::Float64)
    #dz: function that gives derivative of z
    #zi: current state of system before updated rk4Step
    #ti: current time of system before updated rk4Step
    #dt: fixed time step

    k1 = dz(ti, zi)
    k2 = dz(ti + dt/2, zi + dt*k1/2)
    k3 = dz(ti + dt/2, zi + dt*k2/2)
    k4 = dz(ti + dt, zi + dt*k3)
    return zi + 1/6 * (k1 + 2*k2 + 2*k3 + k4)*dt

end

#integrator main function 
function rk4(dz::Function, tarray::Vector{Float64}, z0::Vector{Float64})
    #dz: dericative function
    #tarray: array of times to integrate over (one step per time)
    #z0: initial condition/state
    #returns: z (matrix of states [length(tarray x length(z0))])

    z = zeros(length(tarray), length(z0))
    z[1,:] = z0

    for i = 2:size(z)[1]
        dt = tarray[i] - tarray[i-1]
        z[i,:] = rk4Step(dz, tarray[i-1], z[i-1,:], dt)
    end

    return z
end

##########aero forces##########

#drag force
function getDrag(vRAI_I::Vector{Float64}, Cd::Float64, A::Float64, ρ::Float64)
    #vRAI_I: velocity of rocket wrt to air (m/s)
    #Cd: coefficient of drag
    #A: area normal to velcoty vector
    #ρ: air density (kg/m^3)

    if(norm(vRAI_I) < 1e-15)
        return zeros(3)
    end

    direction = -vRAI_I/norm(vRAI_I);
    mag = .5 * Cd  * ρ * norm(vRAI_I)^2 * A

    return direction * mag

end

#lift force
function getLift(vRAI_I::Vector{Float64}, Cl::Float64, A::Float64, ρ::Float64, q::Vector{Float64})
    #vRAI_I: velocity of rocket wrt to air (m/s)
    #Cl: coefficient of lift
    #A: area normal to velcoty vector
    #ρ: air density (kg/m^3)

    d1 = vRAI_I/norm(vRAI_I); #airflow direction

    if any(isnan.(d1))
        return zeros(3)
    end

    b3_I = rotateFrame([0.0;0.0;1.0], quatInv(q)) #dirction of rocket axis in intertial frame components

    #probably impliment householder rotation soon instead, for now subtract off

    unnormed = b3_I - d1
    liftDirection = unnormed/norm(unnormed)
    mag = .5 * Cl  * ρ * norm(vRAI_I)^2 * A

    return liftDirection * mag
end

#mach number
function calcMach(vRAI_I, h)
    #vRAI_I: velocity of body wrt to air
    #returns: mach number 
    return norm(vRAI_I)/getLocalSoS(h)
end

#velocity of rocket wrt to air
function getVRA(vROI_I::Vector{Float64}, vAOI_I::Vector{Float64})
    #vROI_I: inertial velocity of rocket wrt to ground origin in inertial frame components
    #vAOI_I: inertial velocity of air wrt to ground origin in inertial frame components
    #returns: velocity of rocket wrt to air
    return vROI_I - vAOI_I
end

function calcAoA(vRAI_I::Vector{Float64}, q::Vector{Float64})
    #vRAI_I: velocity of rocket wrt to air
    #q: quaternion that represents frame rotation from inertial to body
    #return AoA: Angle of attack --> Angle between rocket velocty wrt to air and rocket rood chord (centerline)

    #find centerline of rocket (b3) in intertial frame (rotate by conjagate of q)
    b3B = [0.0;0.0;1.0]
    b3I = rotateFrame(b3B, quatInv(q))

    if(norm(vRAI_I) < 1e-13)
        return 0.0
    end
    
    argument = dot(vRAI_I, b3I)/norm(vRAI_I)
    if abs(argument) > 1
        argument = sign(argument)
    end

    return acos(argument)
end

function Ig_solidCylinder(m::Float64, h::Float64, R::Float64)

    return m/12 * [(3*R^2+h^2) 0 0; 
                   0 (3*R^2+h^2) 0; 
                   0 0 (6*R^2)]

end

function getWind(t::Float64, h::Float64)
    #t: time since ignition
    #h: height above sea level
    #return: VAOI_I

    return [0.0;0.0;0.0]

end

function getLocalSoS(h::Float64)
    #h: height above sea level (m)
    #return: speed of sound in m/s

    #retrieved from https://www.engineeringtoolbox.com/elevation-speed-sound-air-d_1534.html

    speeds = [344.1
    342.2
    340.3
    339.3
    338.4
    337.4
    336.4
    335.5
    334.5
    333.5
    332.5
    330.6
    328.6
    326.6
    324.6
    322.6
    320.5]

    altitudes = [-1000
    -500
    0
    250
    500
    750
    1000
    1250
    1500
    1750
    2000
    2500
    3000
    3500
    4000
    4500
    5000]

    interp = LinearInterpolation(altitudes, speeds)
    return interp(h)

end

#updates currentMass element of massData
function updateMassState!(t::Float64, massData::StructArray{massElement}, motorData::Matrix{Float64})
    #t: time
    #massData: !modifies! 
    #motorData

    #rocket specific code
    curThrust, curMass = motorThrustMass(t, motorData, massData.initalMass[2])
    massData.currentMass[2] = curMass
    #update I for new mass

    massData.I[2] = massData.Icalc[2](massData.currentMass[2])

end

function totalAeroForceMoment(t::Float64, z::Vector{Float64}, aeroData::aeroCharacterization, massData::StructArray{massElement})
    #t: time since ignition
    #z: system state
    #aeroData: 
    #massData:
    #returns: total aero force in inertial frame and total aero moment in beta frame

    #useful quantities
    vAOI_I = getWind(t, z[3])
    vRAI_I = getVRA(z[4:6], vAOI_I)
    aoa = calcAoA(vRAI_I, z[7:10])
    mach = calcMach(vRAI_I, z[3])
    cd = getCd(aoa, mach, aeroData)
    cl = getCl(aoa, mach, aeroData)
    A = getA(aoa, aeroData)
    ρ = expAtm_r(z[1:3])

    
    drag_I = getDrag(vRAI_I, cd, A, ρ)
    lift_I = getLift(vRAI_I, cl, A, ρ, z[7:10])

    drag_b = rotateFrame(drag_I, z[7:10])
    lift_b = rotateFrame(lift_I, z[7:10])

    #moment calcs --> only aero
    #transform drag and lift to body frame
    drag_b = rotateFrame(drag_I, z[7:10])
    lift_b = rotateFrame(lift_I, z[7:10])
    #calculate getCOM
    rGR_B = getCOM(massData)
    #calculate COP
    rPR_B = getCOP(aoa, mach, aeroData)
    rPG_B = rPR_B - rGR_B
    #aero force acting at COP
    moment_B = cross(rPG_B, drag_b + lift_b)

    return (drag_I + lift_I), moment_B

end

#state derivative function
function stateDerivative!(t::Float64, z::Vector{Float64}, aeroData::aeroCharacterization, massData::StructArray{massElement}, motorData::Matrix{Float64})
    #t: time
    #z: state
    #aeroData: data describing aero features
    #massData: !modifies! array of mass elements
    #motorData: data describing motor thrust
    #returns: derivative of each element of the state

    updateMassState!(t, massData, motorData) #updates massData given the current time and motor information

    #get aero forces/moments
    totalAero_I, aeroMoment_B = totalAeroForceMoment(t, z, aeroData, massData)

    #get Ig_B
    Ig_B = getIg(massData)

    #accel calc
    #calculate mass
    m = sum(massData.currentMass)
    #Thrust #b3 direction 
    thrustMag = motorThrustMass(t, motorData, massData.initalMass[2])[1]
    thrust_B = [0.0;0.0;thrustMag] #assume thrust along rocket axis
    thrust_I = rotateFrame(thrust_B, quatInv(z[7:10]))
    #Gravity
    grav_I = aGrav(z[1:3])
    
    a_I = (thrust_I + totalAero_I)/m + grav_I

    #handling rocket sitting on pad
    if(t < 1.0 && a_I[3] < 0)
        a_I = [0.0, 0.0, 0.0]
    end

    #based on previous state
    v_I = z[4:6] 
    dqB = .5 * quatProd(z[7:10], [z[11:13];0]) #1/2 * q * [IwB; 0]

    #total moment:
    moment_B = aeroMoment_B
    #calc dqB_B
    dwB_B = Ig_B \ (moment_B - cross(z[11:13], Ig_B * z[11:13])) #moments shit

    return [v_I;a_I;dqB;dwB_B]
end

#location of COM
function getCOM(massData::StructArray{massElement})

    #massData: struct of mass information 
    #return B frame position of total center of mass

    weightedSum = zeros(3)
    for element in massData

        weightedSum = weightedSum + element.currentMass * element.comPos

    end

    return weightedSum/sum(massData.currentMass)

end

#parallel axis theorem calcs, NOT FINISHED
function parallelAxis(Ig::Matrix{Float64}, m::Float64, rQG::Vector{Float64})
    #Ig: inertia tensor about component COM
    #m: total mass of component
    #rQG: position of point you want to know wrt to component COM

    return Ig + m * (I * dot(rQG, rQG) - rQG * rQG')

end

#gets inertia tensor of rocket at current center of mass
function getIg(massData::StructArray{massElement})
    #massData: array of structs describing each mass element
    #return: inertia tensor about center of mass

    total = zeros(3,3)
    rlgR_B = getCOM(massData)

    for element in massData
        rGlg_B = element.comPos - rlgR_B
        total = total + parallelAxis(element.I, element.currentMass, rGlg_B)
    end
    return total
end

#gravity function
function aGrav(rRO_I::Vector{Float64})

    #rRO_I: postition of body wrt to point at sea level directly above or below (radially) launch site (m)
    #return: acceleration due to gravity (vector)

    #assumes sphereical earth. Adds height to mean radius of the earth 
    G = 6.67430e-11 #gravidational constant(kg,m,s)
    Me = 5.97219e24 #mass of the earth (kg)
    rOC_I = [0.0;0.0; 6.3781e6] #radius of the earth (m)

    rRC_I = rOC_I + rRO_I;

    return -G * Me / norm(rRC_I)^3 *rRC_I
    return 

end

#apply quat
function rotateFrame(v::Vector{Float64}, q::Vector{Float64})

    #v: vector in initial frame
    #q: [q1, q2, q3, q] quaternion in vector form with imaginary part as the first 3 elements and real part as the last element
    #returns v in new rotated frame coordinates

    quatVec = [v; 0]
    return quatProd(quatInv(q), quatProd(quatVec, q))[1:3]

end

#quaternion product
function quatProd(q1::Vector{Float64}, q2::Vector{Float64})

    #q1: quaternion [q1,q2,q3,q]
    #q2: quaternion [q1,q2,q3,q]
    #returns: quaternion product q1*q2

    vec = q1[4] * q2[1:3] + q2[4] * q1[1:3] + cross(q1[1:3], q2[1:3])
    real = q1[4] * q2[4] - dot(q1[1:3], q2[1:3])

    return [vec; real]

end

#quaternon inverse
function quatInv(q::Vector{Float64})
    #q: quaternion
    #returns: inverse of q

    return [-q[1:3]; q[4]]

end

#add Row to 2D Matrix
function addRow(matrix, index)

    #matrix: 2D matrix to have row added to
    #index: index of new row
    #return: matrix with new row of zeros at index

    return [matrix[1:(index-1), :]; zeros(size(matrix)[2])'; matrix[index:end, :]]

end

#add Col to 2D Matrix
function addCol(matrix, index)

    #matrix: 2D matrix to have coulumns added to
    #index: position of new column of zeros
    #return: matrix with new column of zeros at index

    return [matrix[:, 1:(index-1)] zeros(size(matrix)[1]) matrix[:, index:end]]

end


function changeTimeData(z, tspan0::Vector{Float64}, dtf::Float64)
    #z: state vector
    #tspan0: original time vector (aligns with z)
    #dtf: desired final timeStep
    #return: tspanf time array corresponding to zf the new data. 

    tspanf = collect(LinRange(0.0, tspan0[end], round(Int, tspan0[end]/dtf) + 1))
    zf = zeros(length(tspanf), size(z)[2])

    for j = 1:length(tspanf)

        #time to find z value
        t = tspanf[j]

        #find index where t is valid
        index = 0
        for i = 2:length(tspan0)
            if(t <= tspan0[i])
                index = i - 1
                break
            end
        end

        zf[j, :] = z[index,:] + (z[index+1,:]-z[index,:])./(tspan0[index+1]-tspan0[index]) .* (t - tspan0[index])
    end

    return tspanf, zf
end

#for use when using Plots
function getQuiverPlot(z::Matrix{Float64}, secondDirection::Int)
    #z: input matrix of state vectors
    #secondDimension: second dimension to be plotted (x or y)
    #gives quiverplot with scatter overlay
    
    #max array of rocket pointing vectors
    b3_I = zeros(size(z)[1], 3)
    for i = 1:size(z)[1]
        b3_I[i,:] = rotateFrame([0,0,1.0], quatInv(z[i, 7:10]))
    end

    #shorten array of b3_I vectors (not as many)
    numToPlot = 50
    array = floor.(Int, collect(LinRange(1, size(z)[1], numToPlot)))
    x = zeros(size(array))
    y = zeros(size(array))
    shortb3_I = zeros(size(array)[1], 3)

    shortIndex = 1
    for index in array
        x[shortIndex] = z[index, secondDirection]
        y[shortIndex] = z[index, 3]
        shortb3_I[shortIndex, :] = b3_I[index, :]
        shortIndex = shortIndex + 1
    end


    #plot  all positions
    plot(z[:,secondDirection], z[:,3])
    #overlay dedensified pointing arrows
    quiver!(x, y, quiver=(shortb3_I[:,secondDirection], shortb3_I[:,3]), aspect_ratio = 1)

end

#for use when using Plots
function getAlignmentPlot(t::Vector{Float64}, z::Matrix{Float64})
    #z: input matrix of state vectors
    #plot: alignmnet of rocket pointing and velocity vector

    displayFrom = 3

    #max array of rocket pointing vectors
    b3_I = zeros(size(z)[1], 3)
    for i = 1:size(z)[1]
        b3_I[i,:] = rotateFrame([0,0,1.0], quatInv(z[i, 7:10]))
    end

    alignment = zeros(size(z)[1])

    for i = 1:size(z)[1]
        alignment[i] = dot(b3_I[i, :], z[i, 4:6])/norm(z[i, 4:6])
    end

    plot(t[displayFrom:end], alignment[displayFrom:end], xlabel = "time(s)", ylabel = "v * b3")

end

#for use when using PyPlot
function getQuiverPlot_py(z::Matrix{Float64}, secondDirection::Int)
    #z: input matrix of state vectors
    #secondDimension: second dimension to be plotted (x or y)
    #gives quiverplot with scatter overlay
    
    #max array of rocket pointing vectors
    b3_I = zeros(size(z)[1], 3)
    for i = 1:size(z)[1]
        b3_I[i,:] = rotateFrame([0,0,1.0], quatInv(z[i, 7:10]))
    end

    #shorten array of b3_I vectors (not as many)
    numToPlot = 30
    array = floor.(Int, collect(LinRange(1, size(z)[1], numToPlot)))
    x = zeros(size(array))
    y = zeros(size(array))
    shortb3_I = zeros(size(array)[1], 3)

    shortIndex = 1
    for index in array
        x[shortIndex] = z[index, secondDirection]
        y[shortIndex] = z[index, 3]
        shortb3_I[shortIndex, :] = b3_I[index, :]
        shortIndex = shortIndex + 1
    end


    #plot  all positions
    plot(z[:,secondDirection], z[:,3])
    #overlay dedensified pointing arrows
    quiver(x, y, shortb3_I[:,secondDirection], shortb3_I[:,3], width=.006)
    axis("equal")
    xlabel("Cross Range (m)")
    ylabel("Local Vertical (m)")

end

function get3DQuiverPlot_py(z::Matrix{Float64})
    #z: input matrix of state vectors
    #gives quiverplot with scatter overlay
    
    #max array of rocket pointing vectors
    b3_I = zeros(size(z)[1], 3)
    for i = 1:size(z)[1]
        b3_I[i,:] = rotateFrame([0,0,1.0], quatInv(z[i, 7:10]))
    end

    #shorten array of b3_I vectors (not as many)
    numToPlot = 30
    array = floor.(Int, collect(LinRange(1, size(z)[1], numToPlot)))
    x = zeros(size(array))
    y = zeros(size(array))
    k = zeros(size(array))
    shortb3_I = zeros(size(array)[1], 3)

    shortIndex = 1
    for index in array
        x[shortIndex] = z[index, 1]
        y[shortIndex] = z[index, 2]
        k[shortIndex] = z[index, 3]
        shortb3_I[shortIndex, :] = b3_I[index, :]
        shortIndex = shortIndex + 1
    end


    #plot  all positions
    plot3D(z[:,1], z[:,2], z[:,3])
    #overlay dedensified pointing arrows
    quiver(x, y, k, shortb3_I[:,1], shortb3_I[:,2], shortb3_I[:,3], length = 10, colors = "red")
    autoscale()
    xlabel("e1 (m)")
    ylabel("e2 (m)")
    zlabel("e3 (m)")

end

#for use when using PyPlot
function getAlignmentPlot_py(t::Vector{Float64}, z::Matrix{Float64})
    #z: input matrix of state vectors
    #plot: alignmnet of rocket pointing and velocity vector

    displayFrom = 3

    #max array of rocket pointing vectors
    b3_I = zeros(size(z)[1], 3)
    for i = 1:size(z)[1]
        b3_I[i,:] = rotateFrame([0,0,1.0], quatInv(z[i, 7:10]))
    end

    alignment = zeros(size(z)[1])

    for i = 1:size(z)[1]
        alignment[i] = dot(b3_I[i, :], z[i, 4:6])/norm(z[i, 4:6])
    end

    plt = plot(t[displayFrom:end], alignment[displayFrom:end])
    ylabel("v * b3")
    xlabel("time(s)")

    return plt

end

function align(tspan, z, offset)
    #tspan: time vector
    #z: state vector as aligns with tspan
    #offset: number of data points to cut off the beginning of z

    zf = z[offset:end, :]
    zf[:,3] = zf[:,3] .- zf[1,3] #reset z3 = 0

    tspanf = tspan[1:end-offset+1]

    return tspanf, zf
end

function penalty(v1, v2)
    #v1: nx1 vector
    #v2: mx1 vector
    #returns sum of the squares of the differences between v1i - v2i for the first min(m,n) components

    length = min(size(v1)[1], size(v2)[1])

    return sum((v1[1:length] .- v2[1:length]).^2)
end

function run_penalty(z0, tspan, aeroData, massData, motorData, flightData)
    #run simulation with given parameters (z0, tspan, aeroData, massData, motorData) and calculate penal

    dz(t, zi) = stateDerivative!(t, zi, aeroData, massData, motorData)
    z = rk4(dz, tspan, z0)

    dt_data = flightData[2,1]
    tspanf, zf = changeTimeData(z, tspan, dt_data)

    tspan_aligned, zf_aligned = align(tspanf, zf, 7)  #hardcoded alignment offset

    return penalty(zf_aligned[:,3], flightData[:,2])

end

function run_plotz3(z0, tspan, aeroData, massData, motorData, flightData)
    #runs simulation based on z0, tspan, aeroData, massData, motorData. 
    #Plots flight data as well

    dz(t, zi) = stateDerivative!(t, zi, aeroData, massData, motorData)
    z = rk4(dz, tspan, z0)

    dt_data = flightData[2,1]
    tspanf, zf = changeTimeData(z, tspan, dt_data)

    tspan_aligned, zf_aligned = align(tspanf, zf, 7)  #hardcoded alignment offset

    pygui(true)
    plot(tspan_aligned, zf_aligned[:,3])
    plot(flightData[1:length(tspanf)+50,1], flightData[1:length(tspanf)+50,2])

end

function run_plotz3(z0, tspan, aeroData, massData, motorData)
    #runs simulation and plots height as a function of time

    dz(t, zi) = stateDerivative!(t, zi, aeroData, massData, motorData)
    z = rk4(dz, tspan, z0)

    pygui(true)
    plot(tspan, z[:,3])

end



function aeroData_Cd_Mach(mach, Cd_Mach)
    #mach: array of mach numbers
    #Cd_mach: corresponding Cd for each mach number in Mach
    #returns: aeroDataSet

    #makes aeroDataSet given coefficients of drag for each mach number. Assumes Cd is constant for each AoA. Only varies with mach number as given 
    #quantities that vary with AoA including A and Cl are hardcoded. 

    aeroData = aeroCharacterization()
    for i = 1:length(mach)
        addData!(aeroData, aeroDataPoint(0.0, mach[i], Cd_Mach[i], 0, -2.8, .02284))
        addData!(aeroData, aeroDataPoint(pi/2, mach[i], Cd_Mach[i], 0, -2.8, .55))
        addData!(aeroData, aeroDataPoint(pi, mach[i], Cd_Mach[i], 0, -2.8, .02284))
    end

    return aeroData
end

#code body
begin

    ######### conventions #################
    #point R is defined at tip of nose cone. All positions relative to rocket are relative to the nosecone. 
    #########parameter definitions###########

    ##Generic Set up##

    #motor info
    thrustCurveFileName = "/Users/Sam/Desktop/Code/6dof/Cesaroni_13628N5600-P.eng"
    motorData = readMotorData(thrustCurveFileName) #s, N
    # initalPropMass = 6.363 #kg

    #read in flight data to compare to 
    flightDataFilePath = "/Users/Sam/Desktop/Code/6dof/SP22CompData.csv"
    flightData = readRRC3Data(flightDataFilePath)

    #dynamic mass properties
    mainBodyIg(m) = Ig_solidCylinder(m, 3.6, .075)
    motorIg(m) = Ig_solidCylinder(m, 1.0, .05)
    massData = StructArray([massElement([0;0;-2.2], 41.036, 41.036, mainBodyIg(41.036), mainBodyIg), massElement([0;0;-3.48], 6.363, 6.363, [1 0 0; 0 1 0; 0 0 1], motorIg)])

    #aero properties (fixed)
    dataSet = aeroCharacterization()
    addData!(dataSet, aeroDataPoint(0.0, 0.0, 0.36, 0, -2.8, .02284))
    addData!(dataSet, aeroDataPoint(pi/2, 0.0, 0.36, .1, -2.8, .55))
    addData!(dataSet, aeroDataPoint(0.0, 2.0, 0.36, 0, -2.9, .02284))
    addData!(dataSet, aeroDataPoint(pi/2, 2.0, 0.36, .1, -2.9, .55))
    
    #state derivative function specific to this rocket + conditions
    dz(t, zi) = stateDerivative!(t, zi, dataSet, massData, motorData)

    #test IC + time
    tspan = collect(LinRange(0.0, 24, 2000))
    r0 = [0.0,0.0, 1400.0]
    v0 = [0.0,0.0,0.0]
    n = [0;1;0]
    θ =  .07
    q0 = [sin(θ/2)*n; cos(θ/2)]
    w0 = zeros(3)
    z0 = [r0;v0;q0;w0]
    #z = rk4(dz, tspan, z0) #solve

    ##  ##  ##  ##  ##

    machArray = [0.0, 0.3, 0.6, 0.9]

    # penalty_Cd_Mach(cd_mach) = run_penalty(z0, tspan, aeroData_Cd_Mach(machArray, cd_mach), massData, motorData, flightData)

    # println("Begin Optimize")
    # lower = ones(length(machArray)) * .2
    # upper = ones(length(machArray)) * .8
    # results = optimize(penalty_Cd_Mach, lower, upper, ones(length(machArray))*.4)
    # optim_cd_mach = Optim.minimizer(results)

    # println(optim_cd_mach)

    run_plotz3(z0, tspan, aeroData_Cd_Mach(machArray, ones(length(machArray))* 0.45), massData, motorData, flightData)

    ############ Past Testing ##########

    #println(run_penalty(z0, tspan, dataSet, massData, motorData, flightData))

    # tspanf, zf = changeTimeData(z, tspan, .05)
    # zf[:,3] = zf[:,3] .- zf[1,3]

    # tspan_aligned, zf_aligned = align(tspanf, zf, 7) #manually aligned by visually trying to match curves at beginning

    # pygui(true)
    # plot(tspan_aligned, zf_aligned[:,3])
    # plot(flightData[:,1], flightData[:,2])

    # p1 = getQuiverPlot(z,2)
    # p2 = getAlignmentPlot(tspan, z)
    # # p3 = plot3d(z[:,1], z[:,2], z[:,3])
    #plot(p1, p2, layout = (1,2))


    # pygui(true)
    # plt = getAlignmentPlot_py(tspan, z)

    # getQuiverPlot_py(z, 2)


    #axis("equal")

    #get3DQuiverPlot_py(z)
    #plot3D(z[:,1], z[:,2], z[:,3])


   
    #test state
    # n = [0;1;0]
    # θ = pi/24 #7.5deg
    # q = [sin(θ/2)*n; cos(θ/2)]
    # ri = [1.0;1.0;100.0]
    # vi = [10.0;0.0;100.0]
    # wi = [.05;0;0]
    # t = 1.5
    # zi = [ri;vi;q;wi]
    
    # println(massData)
    # derivativez = dz(t, zi)
    # println(derivativez)
    # println(massData)

    #testing rotate frame function
    # n = [1;0;0]
    # θ = pi/6
    # q = [sin(θ/2)*n; cos(θ/2)]
    # x3 = [0;0;1]
    # println(rotateFrame(x3,q))

    #testing motor fucnciton 
    #println(motorData)
    # thrust, mass = motorThrustMass(2.3, motorData, initalPropMass)
    # println(thrust)
    # println(mass)

    #testing expAtm function 
    # print("Testing altitude: ")
    # println(expAtm(2))
    # println(expAtm(6))

   

    # println(calcAoA([0;0;100], [0;100;0], q))

    # dataSet = aeroCharacterization()
    # addingSet = aeroDataPoint(.1,.1, 1.0,1.0,1.0,23.0)
    # addData!(dataSet, addingSet)
    # addingSet = aeroDataPoint(.2,.1, 2.0,2.0,2.0,40.0)
    # addData!(dataSet, addingSet)
    # addingSet = aeroDataPoint(.1,.2, 3.0,3.0,3.0,23.0)
    # addData!(dataSet, addingSet)
    # addingSet = aeroDataPoint(.2,.2, 4.0,4.0,4.0,40.0)
    # addData!(dataSet, addingSet)
    # println(dataSet.AoA)
    # println(dataSet.Mach)
    # println(dataSet.Cd)
    # println(dataSet.A)

    # println(getA(.15, dataSet))
    # println(getCd(.19,.15,dataSet))
    # println(getCOP(.15, .12, dataSet))

    

    # println(dataSet.A)
    # println(dataSet.Cd)

    # n = [1;0;0]
    # θ = .05
    # q = [sin(θ/2)*n; cos(θ/2)]

    # vROI_I = [0.0;0.0; 100.0]
    # vAOI_I = [0.0;5.0;0.0]
    # vRAI_I = getVRA(vROI_I, vAOI_I)
    # println(vRAI_I)
    # aoa = calcAoA(vRAI_I, q)
    # println(aoa)
    # println(getDrag(vRAI_I, getCd(aoa, .1, dataSet), getA(aoa, dataSet), expAtm(100)))

    # testingFunction(t,x) = [x[2]; -x[1]]
    # z0 = [1.0; 0.0]
    # tspan = convert(Vector, LinRange(0.0, 10.0, 100))


    # result = rk4(testingFunction, tspan, z0)
 
    # plot(tspan, result)


end
