#include "massElement.cpp"
#include "motor.cpp"
#include <Eigen/Geometry>
#include <vector>
#include "aeroChar_or.cpp"

#define STATE_SIZE 13

typedef Eigen::Matrix<double,STATE_SIZE,1> stateVec;

class rocket{

    private:
        double staticTotalMass;
        double totalMass;
        Eigen::Matrix3d staticIg;
        Eigen::Matrix3d Ig_B;
        Eigen::Vector3d static_rGR_B;
        Eigen::Vector3d rGR_B;
        motor motorData;
        bool burnout;
        aeroChar_or aeroData;


 
    public:

        rocket(){
            staticTotalMass = 0;
            totalMass = 0;
            staticIg = Eigen::Matrix3d::Zero();
            Ig_B = Eigen::Matrix3d::Zero();
            static_rGR_B = Eigen::Vector3d::Zero();
            rGR_B = Eigen::Vector3d::Zero();
            motorData = motor();
            burnout = false;
            aeroData = aeroChar_or();
        }

        //constructors
        rocket(std::vector<massElement> staticMasses, motor m){
            burnout = false;


            //set bulk static mass properties based on staticMasses vector
            motorData = m;
            totalMass = 0;
            staticTotalMass = 0;
            Eigen::Vector3d weightedPos(0,0,0);
            for(int i = 0; i < staticMasses.size(); i++)
            {
                staticTotalMass += staticMasses[i].getMass();
                weightedPos += staticMasses[i].getLocation() * staticMasses[i].getMass();
            }
            static_rGR_B = weightedPos/staticTotalMass;

            staticIg = Eigen::Matrix3d::Zero();
            for(int i = 0; i < staticMasses.size(); i++){

                Eigen::Vector3d riGG_B = staticMasses[i].getLocation() - static_rGR_B; //position of mass wrt rocket static center of mass
                staticIg += parallelAxis(staticMasses[i].getIg(), staticMasses[i].getMass(), staticMasses[i].getLocation());
            }

            //set dynamic mass properties based on motorData, Assume motor is fully loaded
            totalMass = staticTotalMass + motorData.getMass();
            rGR_B = (static_rGR_B * staticTotalMass + motorData.getMass() * motorData.getPos())/totalMass;
            Ig_B = parallelAxis(staticIg, staticTotalMass, rGR_B - static_rGR_B) + parallelAxis(motorData.getIg(), motorData.getMass(), rGR_B - motorData.getPos());

        }

        //functions
        void updateMassState(double time){
            if(time < motorData.getTimes().back()){
                motorData.updateMass(time); //update mass characteristics of motor;
                totalMass = staticTotalMass + motorData.getMass();
                rGR_B = (static_rGR_B * staticTotalMass + motorData.getMass() * motorData.getPos())/totalMass;
                Ig_B = parallelAxis(staticIg, staticTotalMass, rGR_B - static_rGR_B) + parallelAxis(motorData.getIg(), motorData.getMass(), rGR_B - motorData.getPos());
            }
            else if (!burnout){
                totalMass = staticTotalMass;
                rGR_B = static_rGR_B;
                Ig_B = staticIg;
                burnout = true;
            }
        }

        Eigen::Vector3d getCOM(){

            return rGR_B;
        }

        Eigen::Matrix3d getIg(){
            return Ig_B;
        }

        double getMass(){
            return totalMass;
        }

        Eigen::Vector3d getAeroForce(double t, stateVec z){
            return aeroData.getTotalForce_I(t,z);
        }

        Eigen::Vector3d getAeroMoment(double t, stateVec z){

            // returns aeroMoment_B: aeromoments on the rocket written in the B frame

            Eigen::Vector3d rPR_B = aeroData.getCOP(t,z);
            Eigen::Vector3d fP_B = aeroData.getTotalForce_B(t,z);

            return (rPR_B - rGR_B).cross(fP_B);
        }

        Eigen::Vector3d getCOP(double t, stateVec z){
            return aeroData.getCOP(t, z);
        }

        double getThrust(double t){
            return motorData.calcThrust(t);
        }


};

