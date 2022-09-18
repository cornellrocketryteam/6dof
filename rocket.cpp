#include "massElement.cpp"
#include "motor.cpp"
#include <Eigen/Geometry>
#include <vector>

class rocket{

    private:
        std::vector<massElement> staticMasses;
        massElement propellantMass;
        double totalMass;
        motor motorData;
 
    public:

        //constructors
        rocket(std::vector<massElement> masses, motor m){
            motorData = m;
            staticMasses = masses;
            totalMass = 0;
            for(int i = 0; i < std::size(staticMasses); i++)
            {
                totalMass += staticMasses[i].getMass();
            }
            propellantMass = massElement()
        }

        //functions
        void updateMassState(double time){
            totalMass = 0;
            for(int i = 0; i < std::size(staticMasses); i++)
            {
                staticMasses[i].updateMass(time);
                totalMass += staticMasses[i].getMass();
            }
        }

        Eigen::Vector3d getCOM(){

            Eigen::Vector3d weightedPos(0,0,0);

            for(int i = 0; i < std::size(staticMasses); i++)
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

    motor testMotor = motor(thrusts, times, 6.373);

    Eigen::Vector3d propPos (0,0,-3.8);
    massElement propellant = massElement(propPos, &testMotor.calcMass, testMotor.calcMass(0));
}