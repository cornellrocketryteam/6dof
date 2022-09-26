#include <Eigen/Geometry>
#include<iostream>

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
            Ig = Eigen::Matrix3d::Zero();
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

        Eigen::Matrix3d getIg()
        {
            return Ig;
        }


    
};

Eigen::Matrix3d parallelAxis(Eigen::Matrix3d Ig, double m, Eigen::Vector3d rQG){
    //Ig: inertia tensor about component COM
    //m: total mass of component
    //rQG: position of point you want to know wrt to component COM

    return Ig + m * (Eigen::Matrix3d::Identity() * rQG.dot(rQG) - rQG * rQG.transpose());
}

