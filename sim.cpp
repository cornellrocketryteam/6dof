#include "rocket.cpp"
#include "quat.cpp"
#include <fstream>


#define STATE_SIZE 13
#define NUM_STORE 50
#define R_EARTH 6378100
#define G 6.67430e-11
#define M_EARTH 5.97219e24

//set startup-with-shell off

typedef Eigen::Matrix<double,STATE_SIZE,1> stateVec;

void log(std::ofstream* f, double* zlog, double* tspanlog);

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

            //all the meat
            return -z;

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
    Eigen::Vector3d rRC_I (rRO_I[1], rRO_I[2], rRO_I[3] + R_EARTH); 

    return -1 * G * M_EARTH / rRC_I.norm() * rRC_I;
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
    z0 << 1,0,1000,0,0,0,0,0,0,1,0,0,0;

    sim s = sim(testRocket, z0, pow(10,-3), 5000, 25, -70, "export_data.txt");

    s.run();


}

        













