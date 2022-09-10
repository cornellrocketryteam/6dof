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


