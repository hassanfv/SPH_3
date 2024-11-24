#include <cmath>
#include "nr.h"
using namespace std;

void NR::rkqs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &x, const DP htry,
	const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,
	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
	const DP SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
	int i;
	DP errmax,h,htemp,xnew;

	int n=y.size();
	h=htry;
	Vec_DP yerr(n),ytemp(n);
	for (;;) {
		rkck(y,dydx,x,h,ytemp,yerr,derivs); // yerr will be generated. yerr is used to evaluate whether h is a good step size! ytemp will contain evolved y!
		errmax=0.0;
		for (i=0;i<n;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h)); // We do this line to make sure that htemp is never below 0.1 times the previous h (i.e. current h).
		xnew=x+h; // I think xnew is only used to check whether "xnew = x" (meaning h is 0) so that we raise the "stepsize underflow" error !
		if (xnew == x) nrerror("stepsize underflow in rkqs"); // This is used to alert us if h = 0.0 !
	} // When we exit this loop we have the "y" which is evolved for one step h with the best h (i.e. with a h that gives errmax <= 1.0)!
	if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW); // ERRCON: threshold for deciding if step size can be increased significantly.
	else hnext=5.0*h;
	x += (hdid=h); // we add h to x and then set hdid to h !
	for (i=0;i<n;i++) y[i]=ytemp[i]; // We replace y with ytemp which contains the y values with the best h (i.e. h with errmax <= 1.0); see "break" location in code.
}
