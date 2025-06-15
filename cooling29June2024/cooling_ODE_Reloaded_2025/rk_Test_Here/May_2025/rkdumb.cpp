#include "nr.h"

extern Vec_DP *xx_p;
extern Mat_DP *y_p;


// DEFINE y_p as a 2D array !!!!!!!!

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
    y_p[i][0] = v[i];
  }
  xx_p[0] = x1;
  
  x = x1;
  h = (x2 - x1) / nstep;
  
  for (k = 0; k < nstep; k++)
  {
    derivs(x, v, dv);
    rk4_h(v, dv, h, vout, nvar, derivs);
    
    if (x + h == x)
      //nrerror("Step size too small in routine rkdumb");
      cout << "Step size too small in routine rkdumb" << endl;
    
    x += h;
    xx_p[k + 1] = x;
    for (i = 0; i < nvar; i++)
    {
      v[i] = vout[i];
      y_p[i][k+1] = v[i];
    }
  }
    

}

