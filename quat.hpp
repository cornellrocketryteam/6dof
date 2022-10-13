#include <Eigen/Geometry>

Eigen::Quaterniond quatInv(Eigen::Quaterniond q);
Eigen::Vector4d quatInv(Eigen::Vector4d q);
Eigen::Quaterniond quatProd(Eigen::Quaterniond q1, Eigen::Quaterniond q2);
Eigen::Vector4d quatProd(Eigen::Vector4d q1, Eigen::Vector4d q2);
Eigen::Vector4d quatToVec(Eigen::Quaterniond q);
Eigen::Vector3d rotateFrame(Eigen::Vector3d v, Eigen::Quaterniond q);
Eigen::Vector3d rotateFrame(Eigen::Vector3d v, Eigen::Vector4d q);