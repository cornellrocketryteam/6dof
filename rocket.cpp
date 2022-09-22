#include "massElement.cpp"
#include "motor.cpp"
#include <Eigen/Geometry>
#include <vector>
#include "aeroChar_or.cpp"

class rocket{

    private:
        double staticTotalMass;
        double totalMass;
        Eigen::Matrix3d staticIg;
        Eigen::Matrix3d Ig;
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
            Ig = Eigen::Matrix3d::Zero();
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
            Ig = parallelAxis(staticIg, staticTotalMass, rGR_B - static_rGR_B) + parallelAxis(motorData.getIg(), motorData.getMass(), rGR_B - motorData.getPos());

        }

        //functions
        void updateMassState(double time){
            if(time < motorData.getTimes().back()){
                motorData.updateMass(time); //update mass characteristics of motor;
                totalMass = staticTotalMass + motorData.getMass();
                rGR_B = (static_rGR_B * staticTotalMass + motorData.getMass() * motorData.getPos())/totalMass;
                Ig = parallelAxis(staticIg, staticTotalMass, rGR_B - static_rGR_B) + parallelAxis(motorData.getIg(), motorData.getMass(), rGR_B - motorData.getPos());
            }
            else if (!burnout){
                totalMass = staticTotalMass;
                rGR_B = static_rGR_B;
                Ig = staticIg;
                burnout = true;
            }
        }

        Eigen::Vector3d getCOM(){

            return rGR_B;
        }

        Eigen::Matrix3d getIg(){
            return Ig;
        }

        double getMass(){
            return totalMass;
        }

        Eigen::Vector3d getAeroForce(double t, double z[13]){
            return aeroData.getTotalForce(t, z);
        }

        Eigen::Vector3d getCOP(double t, double z[13]){
            return aeroData.getCOP(t, z);
        }




};



int main(){

    double tempArray[] = {6750.0, 5000.0, 4500.0, 3500.0, 2700.0, 0.0};
    std::vector<double> thrusts (std::begin(tempArray), std::end(tempArray));
    double tempArray1[] = {0.019, .5, 1.0, 1.5, 2.0, 2.5};
    std::vector<double> times (std::begin(tempArray1), std::end(tempArray1));
    Eigen::Vector3d propPos (0,0,-3.8);

    motor testMotor = motor(thrusts, times, 6.373, propPos, 1.0, .1);

    double mass = 41.0;
    Eigen::Matrix3d staticIG;
    Eigen::Vector3d staticPos (0,0,-2.8);
    staticIG << mass * Eigen::Matrix3d::Identity();
    massElement staticMass = massElement(staticPos, mass, staticIG);
    std::vector<massElement> masses (1, staticMass);

    rocket testRocket = rocket(masses, testMotor);

    std::cout<<testRocket.getCOM()<<std::endl;
    std::cout<<testRocket.getMass()<<std::endl;
    std::cout<<std::endl;
    testRocket.updateMassState(1.5);
    std::cout<<testRocket.getCOM()<<std::endl;
    std::cout<<testRocket.getMass()<<std::endl;

    
}