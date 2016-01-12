#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

void Hamiltonian(double& H, double* p, double* q);

int main(){
    
    const double dt = 0.0005;
    const double tend = 20*M_PI;
    
    const int N = int(tend/dt);
    const double e = 0.6;
    double p[2], q[2], H, A;
    //Hamiltonian(H, p, q);
    
    ofstream out("Kepler.txt");
    p[0] = 0;
    p[1] = sqrt((1+e)/(1-e));
    q[0] = 1-e;
    q[1] = 0;
    
    
    for (int i = 0; i < N-1; i++){
        
        A = q[0] * q[0] + q[1] * q[1];
        
        p[0] -= dt * q[0] * (pow(A,-3./2));
        p[1] -= dt * q[1] * (pow(A,-3./2));
        q[0] += dt * p[0];
        q[1] += dt * p[1];
        
        Hamiltonian(H, p, q);
        
        out<< i*dt << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
    }
  
  out.close(); 
  return 0;
}    

void Hamiltonian(double& H, double* p, double* q){
    H = 0.5 * (p[0] * p[0] + p[1] * p[1]) - 1. / (sqrt(q[0] * q[0] + q[1] * q[1]));
} 
