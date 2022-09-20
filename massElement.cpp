#include <Eigen/Geometry>

class massElement{

    private:
        Eigen::Vector3d rGR_B;
        double mass;
        Eigen::Matrix3d Ig;

    public:
        //constructor
        massElement(Eigen::Vector3d pos, double m, Eigen::Matrix3d I){
            rGR_B = pos;
            mass = m;
            Ig = I;
        }

        massElement(){
            rGR_B = Eigen::Vector3d (0.0,0.0,0.0);
            mass = 0;
            Ig = Eigen::Matrix3d ();
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
