#include <Eigen/Geometry>

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

        massElement(){
            rGR_B = Eigen::Vector3d (0.0,0.0,0.0);
            mass = 0;
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
