#include "rocket.cpp"
#include "quat.cpp"
#include <fstream>


//constatnts to use
#define STATE_SIZE 13
#define NUM_STORE 100
#define R_EARTH 6378100
#define G 6.67430e-11
#define M_EARTH 5.97219e24
#define W_EARTH 7.292115e-5

//set startup-with-shell off -- Command to use when debugging
//Eigen quick reference guide: https://eigen.tuxfamily.org/dox/group__QuickRefPage.html

typedef Eigen::Matrix<double,STATE_SIZE,1> stateVec;

//function declariations
void log(std::ofstream* f, double* zlog, double* tspanlog);
Eigen::Vector3d aGrav(Eigen::Vector3d rRO_I);
Eigen::Vector3d aCorriolis(double vROI_I_3, double lat);
Eigen::Vector3d aCentrifugal(double lat);

class sim{

    private:

        rocket lv;
        stateVec z0;  //initial conditions
        double dt;
        double t0;
        int n;
        stateVec zout; // all states to return
        std::string fileName;
        double latitude;
        double longitude;

    public:

        sim(rocket r, stateVec zinit, double deltt, int num, double lat, double longi, std::string file){

            lv = r;
            dt = deltt;
            n = num;
            z0 = zinit;
            latitude = lat;
            longitude = longi;
            t0 = 0;
            fileName = file;

        }

        stateVec stateDerivative(double t, stateVec z){

            //compeonents of stateDerivative to build
            Eigen::Vector3d a_I;
            Eigen::Vector4d dqB (0,0,0,0);
            Eigen::Vector3d dwB_B (0,0,0);

            //aero effects
            Eigen::Vector3d fP_I = lv.getAeroForce(t,z); //total aeroforces in the body frame
            Eigen::Vector3d mP_B = lv.getAeroMoment(t,z);

            a_I << aGrav(z.head<3>());
            a_I[2] = a_I[2] + (lv.getThrust(t)/lv.getMass());
            a_I = a_I + (lv.getAeroForce(t,z))/lv.getMass() + aCorriolis(z[5], latitude) + aCentrifugal(latitude);

            if(t < 0.1 && a_I[2] < 0){
                a_I << 0,0,0;
            }
            
            Eigen::Vector4d w_B_quat;
            w_B_quat << z.segment<3>(10), 0;
            dqB = 0.5 * quatProd(z.segment<4>(6), w_B_quat);
            dwB_B = lv.getIg().householderQr().solve(mP_B - z.segment<3>(10).cross(lv.getIg() * z.segment<3>(10)));

            //assembling components to return
            stateVec ret;
            ret << z.segment<3>(3), a_I, dqB, dwB_B;
            return ret;

        }

        stateVec rk4Step(stateVec z_vec, double ti, double dt){

            stateVec k1, k2, k3, k4;
            k1 = this->stateDerivative(ti, z_vec);
            k2 = this->stateDerivative(ti + 0.5 * dt, z_vec + 0.5 * dt * k1);
            k3 = this->stateDerivative(ti + 0.5 * dt, z_vec + 0.5 * dt * k2);
            k4 = this->stateDerivative(ti + dt, z_vec + dt * k3);

            return z_vec + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4) * dt;
        }

        void run(){

            //arrays to temp store
            double* zlog = new double[STATE_SIZE * NUM_STORE]; //z0
            double* tspanlog = new double[NUM_STORE];
    
            stateVec current = z0;
            double currentTime = t0;
            std::ofstream file;
            file.open(fileName, std::ofstream::out | std::ofstream::trunc);

            //put first state in the storeArray
            tspanlog[0] = currentTime;
            for(int i = 0; i < STATE_SIZE; i++){
                zlog[i] = z0[i];
            }

            int j = 1; //counter for states in array

            for(int i = 0; i < n; i++){

                lv.updateMassState(currentTime);

                current = rk4Step(current, currentTime, dt);
                currentTime += dt;

                //fill zlog with current state
                tspanlog[j] = currentTime;
                for(int k = 0; k < STATE_SIZE; k++){
                    zlog[j*STATE_SIZE + k] = current[k];
                }

                //incriment log array
                j++;

                if (j == NUM_STORE || i == (n-1))
                {
                    j = 0;
                    //write array to file
                    log(&file, zlog, tspanlog);
                }
            }
            file.close();
            zout = current;
        }

        stateVec run(stateVec z0_push, double t0_push){
            z0 = z0_push;
            t0 = t0_push;
            run();
            return zout;
        }

        //get & set functions
        void set_dt(double dt_set){
            dt = dt_set;
        }

        stateVec get_zout(){
            return zout;
        }
};

void log(std::ofstream* f, double* zlog, double* tspanlog){

    for(int k = 0; k < NUM_STORE; k++){
        *f << tspanlog[k] << ",";
        tspanlog[k] = 0;
        for(int l = 0; l < STATE_SIZE; l++)
        {
            *f << zlog[(k*STATE_SIZE) + l] << ",";
            zlog[(k*STATE_SIZE) + l] = 0.0;
        }
        *f << std::endl;
    }

}

Eigen::Vector3d aGrav(Eigen::Vector3d rRO_I){
    //rRO_I: postition of body wrt to point at sea level directly above or below (radially) launch site (m)
    //return: acceleration due to gravity (vector)

    //assumes sphereical earth. Adds height to mean radius of the earth 
    Eigen::Vector3d rRC_I (rRO_I[0], rRO_I[1], rRO_I[2] + R_EARTH); 

    return -1 * G * M_EARTH / pow(rRC_I.norm(),3) * rRC_I;
}
Eigen::Vector3d aCorriolis(double vROI_I_3, double lat){
    return Eigen::Vector3d(-2 * W_EARTH * vROI_I_3 * cos(lat),0,0);
}

Eigen::Vector3d aCentrifugal(double lat){
    //lat: latitiude of launch site

    //assumes displacement from surface is negligible relative to the R_EARTH
    return Eigen::Vector3d(0, -sin(lat) * pow(W_EARTH,2) * R_EARTH, cos(lat) * pow(W_EARTH,2) * R_EARTH);
}

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

    stateVec z0 (STATE_SIZE);
    z0 << 0,0,1000,0,0,0,0,0,0,1,0,0,0;

    sim s = sim(testRocket, z0, 5.0 * pow(10,-2), 999, 25, -70, "export_data.txt");

    s.run();


}

        













