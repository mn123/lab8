// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double t);
void RKstep(double* const yn, const double* const y0, const double t, const double dt, double* k1, double* k2, double* k3, double* k4, const int dim);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 2;
	double t, p0;
	const double L = 20, dt = 0.15;
for (p0=0; p0<5; p0 +=0.01){
  t=0;
  double y0[dim] = {p0, 0};
	double yn[dim];
  double k1[dim], k2[dim], k3[dim], k4[dim];
  out << t << "\t" << y0[0] << "\t" << y0[1] << endl;
	while(t<=L)
	{
		t += dt;
		RKstep(yn, y0, t, dt,k1,k2,k3,k4,dim);
		if (y0[1]>0 && yn[1]<0)
		  break;
		
		  
                for(int i=0; i<dim; i++) y0[i] = yn[i];
		out << t << "\t" << y0[0] << "\t" << y0[1] << endl;
	}
	
	double lth=0, rth=1, th;
	double b1,b2,b3,b4;
	t -= dt;
	for (int i=0; i<10; i++){
	  th=(lth+rth)/2.0;
	  b1=th-3*th*th/2.0+2*th*th*th/3.0;
	  b2=th*th-2*th*th*th/3.0;
	  b3=b2;
	  b4=-th*th/2.0+2*th*th*th/3.0;
	  //yn[0]=y0[0]+dt*(b1*k1[0]+b2*k2[0]+b3*k3[0]+b4*k4[0]);
	  yn[1]=y0[1]+dt*(b1*k1[1]+b2*k2[1]+b3*k3[1]+b4*k4[1]);
	
	  if (yn[1]>0)
	    lth=th;
	  else rth=th;
	}
	out << t+dt*th << "\t" << yn[0] << "\t" << yn[1] << endl;
	out.close();
	cout << p0 << "\t" << t+dt*th << endl;
  
}
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double t, const double dt, double* k1, double* k2, double* k3, double* k4, const int dim)
{
	for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, t);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dt * k1[i];
	f(k2, t+0.5*dt);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dt * k2[i];
	f(k3, t+0.5*dt);

	for(int i=0;i<dim; i++) k4[i] = y0[i] + dt * k3[i];
	f(k4,  t+dt);

	for(int i=0;i<dim; i++)
	yn[i] = y0[i] + 1./6.*dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double t)
{
	double y[2] = { y0[0], y0[1] };

	y0[0] = y0[1];
	y0[1] = -y[0]/sqrt(1+y[0]*y[0]);

}
