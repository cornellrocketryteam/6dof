#include <Eigen/Dense>
#include <Eigen/Geometry>

Eigen::Quaterniond quatInv(Eigen::Quaterniond q){

    Eigen::Quaterniond returnQuat;
    returnQuat.vec() = -1 * q.vec();
    returnQuat.w() = q.w();

    return returnQuat;
}

Eigen::Quaterniond quatProd(Eigen::Quaterniond q1, Eigen::Quaterniond q2) {
    
    Eigen::Quaterniond result;
    result.setIdentity();

    result.vec() = q1.w() * q2.vec() + q2.w() * q1.vec() + q1.vec().cross(q2.vec());
    result.w() = q1.w() * q2.w() - q1.vec().dot(q2.vec());

    return result;
}

Eigen::Vector3d rotateFrame(Eigen::Vector3d v, Eigen::Quaterniond q){

    Eigen::Quaterniond quatVec;
    quatVec.vec() = v;
    quatVec.w() = 0.0;

    return quatProd(quatInv(q), quatProd(quatVec, q)).vec();

}