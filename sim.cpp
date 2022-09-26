#include "rocket.cpp"
#include "quat.cpp"


const int stateSize = 13;

typedef Eigen::Matrix<double,stateSize,1> stateVec;
stateVec rk4Step(stateVec (*dz)(double, stateVec), stateVec z_vec, double ti, double dt);

class sim{

    private:

        rocket bigRed1;
        std::array<double, stateSize> z0;  //initial conditions
        std::vector<double> t;
        std::vector<double> zout; // all states to return

    public:

        sim(rocket r, std::vector<double> tspan,  std::array<double, stateSize> zinit){
            bigRed1 = r;
            t = tspan;
            zout.reserve(t.size() * stateSize);
            z0 = zinit;
        }

        std::array<double, stateSize> stateDerivative(double t, std::array<double, stateSize>  z){

        }

        


};

std::array<double, stateSize> mult(double c, std::array<double, stateSize> v){
    std::array<double, stateSize> r;
    for(int i = 0; i < stateSize; i++){
        r[i] = c * v[i];
    }
    return r;
}

std::array<double, stateSize> add(std::array<double, stateSize> v1, std::array<double, stateSize> v2){

    std::array<double, stateSize> r;
    for(int i = 0; i < stateSize; i++){
        r[i] = v1[i] + v2[i];
    }

    return r;
}

stateVec rk4Step(stateVec (*dz)(double, stateVec), stateVec z_vec, double ti, double dt){

    stateVec k1, k2, k3, k4;
    k1 = dz(ti, z_vec);
    k2 = dz(ti + 0.5 * dt, z_vec + 0.5 * dt * k1);
    k3 = dz(ti + 0.5 * dt, z_vec + 0.5 * dt * k2);
    k4 = dz(ti + dt, z_vec + dt * k3);

    return z_vec + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4) * dt;
}



