
#include <Eigen/Dense>

#define STATE_SIZE 13

typedef Eigen::Matrix<double,STATE_SIZE,1> stateVec;

class aeroChar_or{

    private:




    public:


        //constructors


        //functions 
        Eigen::Vector3d getCOP(double t, stateVec z){

            //returns rPG_B. Body frame location of center of pressure

            Eigen::Vector3d rPR_B (0,0,0);

            return rPR_B;

        }

        Eigen::Vector3d getTotalForce_B(double t, stateVec z){

            //returns fP_B: total body aerodynamic forces acting at the center of pressure. 

            Eigen::Vector3d fP_B (0,0,0);

            return fP_B;
        }

        Eigen::Vector3d getTotalForce_I(double t, stateVec z){

            //returns fP_I: total body aerodynamic forces acting at the center of pressure given in inertial frame compoennts. 

            Eigen::Vector3d fP_I (0,0,0);

            return fP_I;
        }


};