#include<vector>
#include<iostream>
#include<Eigen/Dense>


class motor{
    private:
        std::vector<double> thrusts;
        std::vector<double> times;
        double initPropMass;
        Eigen::Vector3d rPR_B;


    public:

        //constructor
        motor(std::vector<double> thrustArray, std::vector<double> timeArray, double propMass, Eigen::Vector3d r)
        {
            thrusts = thrustArray;
            times = timeArray;
            initPropMass = propMass;
            rPR_B = r;
        }

        motor() = default;

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

};


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

int main(){

    double tempArray[] = {6750.0, 5000.0, 4500.0, 3500.0, 2700.0, 0.0};
    std::vector<double> thrusts (std::begin(tempArray), std::end(tempArray));
    double tempArray1[] = {0.019, .5, 1.0, 1.5, 2.0, 2.5};
    std::vector<double> times (std::begin(tempArray1), std::end(tempArray1));

    motor testMotor = motor(thrusts, times, 6.373);

    std::cout << testMotor.calcThrust(1.365394852) << std::endl;

    std::cout << testMotor.calcThrust(10.0) << std::endl;
    std::cout << testMotor.calcMass(1.2) << std::endl;
    std::cout << testMotor.calcMass(10.0) << std::endl;
    
}