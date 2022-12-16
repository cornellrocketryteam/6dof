using LinearAlgebra

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

function quatMag(q::Vector{Float64})
    #q: quaternion (or hopefully)
    #returns: magnitude (should always be 1 if a true quaternion)

    return norm(q)

end

function quatError(q1::Vector{Float64}, q2::Vector{Float64})
    #q1: quaternion [q1,q2,q3,q]
    #q2: quaternion [q1,q2,q3,q]
    #returns: the error quaternion q1*q2^-1

    return quatProd(q1, quatInv(q2))

end

function rev2quatError(ds::Vector{Float64}, f::Float64, a::Float64)
    #ds: generalized Rodrigues error vector
    #f: scale factor
    #a: parameter from 0 to 1
    #returns: corresponding error quaternion

    #Crassidis, John. (2006). Sigma-Point Kalman Filtering for Integrated GPS and Inertial Navigation. Aerospace and Electronic Systems, IEEE Transactions on. 42. 750 - 756. 10.1109/TAES.2006.1642588. 

    ds_norm = norm(ds)
    dq4 = (-a * ds_norm^2 + f * sqrt(f^2 + (1-a^2) * ds_norm^2))/(f^2 + ds_norm^2)
    dqvec = (a + dq4)/f * ds

    return [dqvec; dq4]

end

function quatError2rev(dq::Vector{Float64}, f::Float64, a::Float64)

    return (f/(a + dq[4])) * dq[1:3]

end


