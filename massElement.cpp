#include <Eigen/Geometry>
#include "motor.cpp"

class massElement{

    private:
        Eigen::Vector3d rGR_B;
        double mass;

    public:
        //constructor
        massElement(Eigen::Vector3d pos, double m){
            rGR_B = pos;
            mass = m;
        }

        //gets
        double getMass(){
            return mass;
        }
        void setMass(double m){
            mass = m;
        }
        Eigen::Vector3d getLocation(){
            return rGR_B;
        }

    
};

int main(){

    
}
