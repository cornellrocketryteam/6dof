#include<vector>
#include<iostream>
#include<Eigen/Dense>
#include <math.h> 


Eigen::Matrix3d IgSolidCylinder(double m, double h, double R);


class motor{

    private:
        std::vector<double> thrusts;
        std::vector<double> times;
        double initPropMass;
        double mass;
        Eigen::Matrix3d propIg;
        Eigen::Vector3d rPR_B;
        double R, L;

    public:

        //constructor
        motor(std::vector<double> thrustArray, std::vector<double> timeArray, double propMass, Eigen::Vector3d r, double L, double R)
        {
            thrusts = thrustArray;
            times = timeArray;
            initPropMass = propMass;
            mass = initPropMass;
            rPR_B = r;
            propIg = IgSolidCylinder(propMass, L, R);
            R = R;
            L = L;

        }

        motor(){
            thrusts = std::vector<double>();
            times = std::vector<double>();
            initPropMass = 0;
            mass = 0;
            rPR_B = Eigen::Vector3d (0,0,0);
            propIg = Eigen::Matrix3d::Zero();
            R = 0;
            L = 0;
        }

        void updateMass(double time)
        {
            mass = calcMass(time);
            propIg = IgSolidCylinder(mass, L, R);
        }

        //function to calculate mass of propellant at any given time
        double calcMass(double time){
            double tend = times.back();

            if(time >= tend){
                return 0.0;
            }
            else{
                return initPropMass * (1 - time/tend);
            }
        }


        //function to calculate thrust of motor at any given time
        double calcThrust(double time){

            double tend = times.back();
            //std::cout << tend;

            if(time >= tend){
                return 0.0;
            }
            else{
                int i = 0;
                while(i < times.size() - 1 && time > times[i + 1]){
                    i++;
                }
                return thrusts[i] + (thrusts[i+1]-thrusts[i])/(times[i+1]-times[i]) * (time - times[i]);
            }
            
        }

        std::vector<double> getTimes()
        {
            return times;
        }

        std::vector<double> getThrusts()
        {
            return thrusts;
        }
        
        double getInitPropMas()
        {
            return initPropMass;
        }

        double getMass()
        {
            return mass;
        }

        Eigen::Vector3d getPos()
        {
            return rPR_B;
        }

        Eigen::Matrix3d getIg()
        {
            return propIg;
        }

};

Eigen::Matrix3d IgSolidCylinder(double m, double h, double R){
    Eigen::Matrix3d Ig;
    Ig << (3 * pow(R,2) + pow(h,2)),0,0,
           0,(3 * pow(R,2) + pow(h,2)),0,
           0,0,(6* pow(R,2));

    return Ig;
}


double calcMass(motor m, double time){
    std::vector<double> times = m.getTimes();
    std::vector<double> thrusts = m.getThrusts();
    double initPropMass = m.getInitPropMas();

    double tend = times.back();

    if(time >= tend){
        return 0.0;
    }
    else{
        return initPropMass * (1 - time/tend);
    }
}



// int main(){

//     TESTING
//     double tempArray[] = {6750.0, 5000.0, 4500.0, 3500.0, 2700.0, 0.0};
//     std::vector<double> thrusts (std::begin(tempArray), std::end(tempArray));
//     double tempArray1[] = {0.019, .5, 1.0, 1.5, 2.0, 2.5};
//     std::vector<double> times (std::begin(tempArray1), std::end(tempArray1));
//     Eigen::Vector3d r; r << 0,0,-3.0;

//     motor testMotor = motor(thrusts, times, 6.373, r, 1.0, .1);

//     std::cout << testMotor.calcThrust(1.365394852) << std::endl;
//     std::cout << testMotor.calcThrust(10.0) << std::endl;
//     std::cout << testMotor.calcMass(1.2) << std::endl;
//     std::cout << testMotor.calcMass(10.0) << std::endl;
    
// }