#include <Eigen/Dense>
#include <Eigen/Geometry>

Eigen::Quaterniond quatInv(Eigen::Quaterniond q){

    Eigen::Quaterniond returnQuat;
    returnQuat.vec() = -1 * q.vec();
    returnQuat.w() = q.w();

    return returnQuat;
}

Eigen::Vector4d quatInv(Eigen::Vector4d q){
    return Eigen::Vector4d (-q[0], -q[1], -q[2], q[3]);
}

Eigen::Quaterniond quatProd(Eigen::Quaterniond q1, Eigen::Quaterniond q2) {
    
    Eigen::Quaterniond result;
    result.setIdentity();

    result.vec() = q1.w() * q2.vec() + q2.w() * q1.vec() + q1.vec().cross(q2.vec());
    result.w() = q1.w() * q2.w() - q1.vec().dot(q2.vec());

    return result;
}

Eigen::Vector4d quatProd(Eigen::Vector4d q1, Eigen::Vector4d q2){

    Eigen::Vector4d vec;
    Eigen::Vector3d q1v;
    Eigen::Vector3d q2v;
    q1v << q1.head<3>();
    q2v << q2.head<3>();
    vec << q1[3] * q2v + q2[3] * q1v + q1v.cross(q2v), q1[3] * q2[3] - q1v.dot(q2v);

    return vec;
    
}

Eigen::Vector4d quatToVec(Eigen::Quaterniond q){
    return Eigen::Vector4d (q.vec()[0], q.vec()[1], q.vec()[2], q.w());
}

Eigen::Vector3d rotateFrame(Eigen::Vector3d v, Eigen::Quaterniond q){

    Eigen::Quaterniond quatVec;
    quatVec.vec() = v;
    quatVec.w() = 0.0;

    return quatProd(quatInv(q), quatProd(quatVec, q)).vec();

}

Eigen::Vector3d rotateFrame(Eigen::Vector3d v, Eigen::Vector4d q){

    Eigen::Vector4d quatVec;
    quatVec << v, 0;

    return quatProd(quatInv(q), quatProd(quatVec,q)).head<3>();
}

// int main(){

//     Eigen::Vector4d q1(.70710678118,0,0,.70710678118);
//     Eigen::Vector4d q2(0,0.5,0,0.5);
//     Eigen::Vector3d v(0,0,1);


//     std::cout << rotateFrame(v, q1) <<std::endl;
// }