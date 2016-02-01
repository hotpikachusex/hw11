#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* f1, cmplx* f0, const int Nx, const double k, const double dx, const double xmin,
	  const double dt, double x);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40.0;
  const double xmax = 40.0;
	const double Tend = 10.0 * M_PI;
	const double dx =(xmax - xmin) / (Nx-1);
	const double dt = dx / 100.0;
	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);
	double x;
	const double lambda = 10.0;
	const double omega  = 0.2;
	const double k      = omega * omega;
	const double alpha  = pow(k , 1/4.0);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
	
		  step(psi1, psi0, Nx, k, dx, xmin, dt, x);
		  
		  h = psi0;
		  psi0 = psi1;
		  psi1 = h;
		  t+=dt;
		  
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
	
	delete[] psi0;
	delete[] psi1;
	return 0;
}
//-----------------------------------
void step(cmplx* f1, cmplx* f0, const int Nx, const double k, const double dx, const double xmin,
	  const double dt, double x){
  
  cmplx* d 	= new cmplx[Nx];
  cmplx* u 	= new cmplx[Nx];
  cmplx* l 	= new cmplx[Nx];
  cmplx* rechts = new cmplx[Nx];
  double* V 	= new double[Nx];
  
  for(int i=0;i<Nx;i++){ 
  x = xmin + i * dx;
  V[i] = 0.5 * k * x * x;
  d[i] = cmplx(1.0 , (dt / (2.0 * dx * dx) + dt / (2.0) * V[i]));
  u[i] = cmplx(0.0 , -(dt / (4,0 * dx * dx)));
  l[i] = cmplx(0.0 , -(dt / (4,0 * dx * dx)));
  }
  
  //right side A_stern * Psi_0
  //---------------------------------------------------------------------  

  //erster Wert  = d_0_stern * Psi_0[0] + u_0_stern * Psi_0[1]			
  rechts[0] 	 = cmplx(1.0 , (- (dt / (2.0 * dx * dx) + dt / (2.0) * V[0]))) * f0[0] + cmplx(0.0 , dt / (4,0 * dx * dx)) * f0[1];
 
  // alle Werte dazwischen
  for(int i = 1; i < Nx - 1;i++){
    x = xmin + i * dx;
    rechts[i] = cmplx(0.0 , dt / (4,0 * dx * dx)) * f0[i-1] + cmplx(1 , (- (dt / (2.0 * dx * dx) + dt / (2.0) * V[i]))) * f0[i]
	      + cmplx(0.0 , dt / (4,0 * dx * dx)) * f0[i+1];
  }
  
  //letzter Wert
  
  rechts[Nx-1] = cmplx(0.0 , dt / (4,0 * dx * dx)) * f0[Nx-2] + cmplx(1 , (- (dt / (2.0 * dx * dx) + dt / (2.0) * V[Nx-1]))) * f0[Nx-1];
  
  // forward substitution
  //--------------------------------------------------------------------
  
  for(int j = 0; j < Nx-1; j++){
   d[j+1]  	-= u[j] * l[j+1] / d[j];
   rechts[j+1] 	-= rechts[j] * l[j+1] / d[j]; 
   }
   
  //backward substitution
  //--------------------------------------------------------------------
  
  //letzter Wert
  f1[Nx-1] = rechts[Nx-1] / d[Nx-1];
  
  //alle anderen Werte einschließlich des Ersten
  for(int j = Nx-2; j >= 0; j--){
    f1[j] = (rechts[j] - u[j] * f1[j+1]) / d[j];
  }
   
  
 
  delete[] d;
  delete[] u;
  delete[] l;
  delete[] V;
  delete[] rechts;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
 
