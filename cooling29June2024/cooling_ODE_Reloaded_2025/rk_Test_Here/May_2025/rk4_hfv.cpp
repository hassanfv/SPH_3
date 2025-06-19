#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;


const int nvar = 3; // Lorenz chaotic problem !
const int nstep = 1000;

float sigma = 10.0;
float beta = 8.0 / 3.0;
float rho = 28.0;

float *xx_p = new float[nstep+1];
float *y_p = new float[nvar * (nstep+1)]; // y_p[nvar][nstep+1]... nvar is nrow.


//----- derivs
void derivs(const float x, float *y, float *dydx)
{
  dydx[0] = sigma * (y[1] - y[0]);
  dydx[1] = y[0] * (rho - y[2]) - y[1];
  dydx[2] = y[0] * y[1] - beta * y[2];
}


//----- rk4_h
void rk4_h(float *y, float *dydx, const float x, const float h, float *yout, const int N,
           void derivs(const float, float*, float*))
{
  float xh, hh, h6;
  
  float dym[N], dyt[N], yt[N];
  
  hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;
	
	int i;
	
	for (i = 0; i < N; i++) yt[i] = y[i] + hh * dydx[i]; // dydx = k1
	derivs(xh, yt, dyt); // dyt = k2
	
	for (i = 0; i < N; i++) yt[i] = y[i] + hh * dyt[i]; // dyt = k2
	derivs(xh, yt, dym); // dym = k3
	
	for (i = 0; i < N; i++)
	{
		yt[i] = y[i] + h * dym[i]; // dym = k3
		dym[i] += dyt[i]; // k2 + k3. Note that dyt is k2 from above and we are adding it to dym which is k3. Then dym contains both k2 and k3 !
	}
	
	derivs(x + h, yt, dyt); // Now dyt is overwritten as k4 because it was already added to dym as k2 so it is not needed any more and we can asigned it to be k4!
	
	for (i = 0; i < N; i++)
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]); // in paranthesis we have, in order (k1 + k4 + 2.0*(k2+k3))
}


//---- rkdumb_h
void rkdumb_h(float *vstart, const float x1, const float x2, int nvar, int nstep,
              void derivs(const float, float*, float*))
{
  float x, h;
  float *v, *vout, *dv;
  
  v = new float[nvar];
  vout = new float[nvar];
  dv = new float[nvar];
  
  int i, k;
  for (i = 0; i < nvar; i++)
  {
    v[i] = vstart[i];
    y_p[i * (nstep + 1) + 0] = v[i];
  }
  xx_p[0] = x1;
  
  x = x1;
  h = (x2 - x1) / nstep;
  
  // main for loop
  for (k = 0; k < nstep; k++)
  {
    derivs(x, v, dv);
    rk4_h(v, dv, x, h, vout, nvar, derivs);
    
    if (x + h == x)
      //nrerror("Step size too small in routine rkdumb");
      cout << "Step size too small in routine rkdumb" << endl;
    
    x += h;
    xx_p[k + 1] = x;
    for (i = 0; i < nvar; i++)
    {
      v[i] = vout[i];
      y_p[i * (nstep + 1) + (k + 1)] = v[i];
    }
  }
}


int main()
{

  float y0[3] = {-8.0, 8.0, 27.0};
  const float x1 = 0.0;
  const float x2 = 10.0;
  
   auto start = std::chrono::high_resolution_clock::now();
  
  rkdumb_h(y0, x1, x2, nvar, nstep, derivs);
  
  auto end = std::chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = end - start;
  cout << "Elapsed time: " << elapsed.count() << " seconds" << endl;
  
  
  ofstream outfile("out.csv"); // use pplot.py to plot the result !
  
  outfile << "x,y" << endl;
  
  for (int j = 0; j < nstep+1; j++)
  {
    int n_row_1 = 1;
    int n_row_2 = 2;
    outfile << y_p[n_row_1 * (nstep+1) + j] << "," << y_p[n_row_2 * (nstep+1) + j] << endl;
  }
  
  outfile.close();

}


