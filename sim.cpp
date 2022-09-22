#include "rocket.cpp"
#include "quat.cpp"


const int stateSize = 13u;


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

        std::array<double, stateSize> rk4Step(double ti, double dt, std::array<double, stateSize> zi){

            std::array<double, stateSize> k1,k2,k3,k4;
            k1 = this->stateDerivative(ti, zi);
            k2 = this->stateDerivative(ti + dt/2, add(zi, mult(0.5 * dt, k1)));
            k3 = this->stateDerivative(ti + dt/2, add(zi, mult(0.5 * dt, k2)));
            k4 = this->stateDerivative(ti + dt, add(zi, mult(dt, k3)));
            
            return 
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



