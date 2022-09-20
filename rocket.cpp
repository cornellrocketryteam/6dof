#include "massElement.cpp"
#include "motor.cpp"
#include <Eigen/Geometry>
#include <vector>

class rocket{

    private:
        std::vector<massElement> staticMasses;
        massElement propellant;
        double totalMass;
        Eigen::Matrix3d Ig;
        Eigen::Vector3d rGR_B;
        motor motorData;
 
    public:

        //constructors
        rocket(std::vector<massElement> masses, motor m){
            motorData = m;
            staticMasses = masses;
            totalMass = 0;
            for(int i = 0; i < staticMasses.size(); i++)
            {
                totalMass += staticMasses[i].getMass();
            }
            propellant = massElement(m.getPos(), m.getInitPropMas());
            totalMass += propellant.getMass();
        }

        //functions
        void updateMassState(double time){
            totalMass -= propellant.getMass();
            propellant.setMass(motorData.calcMass(time));
            totalMass += propellant.getMass();
        }

        Eigen::Vector3d getCOM(){

            Eigen::Vector3d weightedPos(0,0,0);

            for(int i = 0; i < staticMasses.size(); i++)
            {
               weightedPos += staticMasses[i].getLocation() * staticMasses[i].getMass();
            }

            return weightedPos/totalMass;
        }


};


int main(){

    double tempArray[] = {6750.0, 5000.0, 4500.0, 3500.0, 2700.0, 0.0};
    std::vector<double> thrusts (std::begin(tempArray), std::end(tempArray));
    double tempArray1[] = {0.019, .5, 1.0, 1.5, 2.0, 2.5};
    std::vector<double> times (std::begin(tempArray1), std::end(tempArray1));
    Eigen::Vector3d propPos (0,0,-3.8);

    motor testMotor = motor(thrusts, times, 6.373, propPos);

    
}