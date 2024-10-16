%%writefile test.cu
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <random>
#include <tuple>
#include "hfvCLibs_v7_hfv_v8_gradRho_v4t.h"
#include "bh_tree_iteration_v2.h"
#include "ngb_v5.1.h"
#include <cstdlib> // This is ONLY used for the "exit(0)" function !!

// hcool function updated (10 Jan 2024).
// Adaptive nSplit implemented (5 Jan 2024).
// Added the isothermal gravitational field acceleration. (24 May 2023).
// Added the reading of the params.txt file and updated the IC reading file section and function. (22 May 2023).

using namespace std;


int main()
{

  float dt = 2e-7; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This is only the first time step !!

  const float Nngb_f = 64.0f; // used in smoothing func.
  const int Nngb = 64;
  const int Ndown = Nngb - 5;
  const int Nup = Nngb + 5;
  const float coeff = 0.005f; // used for smoothing length.
  
  //const float kpc_in_cm = 3.086e21;
  
  //************************************************************************
  //********************* Reading CLOUDY Cooling File **********************
  //************************************************************************
  float XH = 0.7;
  float kB = 1.3807e-16;
  float mH = 1.673534e-24;
  float gamma = 5.0f/3.0f;
  
  ifstream binFile("HCoolMu.bin", ios::binary);
  if (!binFile.is_open()) 
  {
    cerr << "Failed to open HCoolMu.bin file" << endl;
    return 1;
  }

  // Scalar values
  int N_nH, N_T, N_r, N_NH;
  float nHLowBound, nHUpBound, TLowBound, TUpBound, rLowBound, rUpBound, NHLowBound, NHUpBound;

  // Read scalar values
  binFile.read(reinterpret_cast<char*>(&N_nH), sizeof(N_nH));
  binFile.read(reinterpret_cast<char*>(&N_T), sizeof(N_T));
  binFile.read(reinterpret_cast<char*>(&N_r), sizeof(N_r));
  binFile.read(reinterpret_cast<char*>(&N_NH), sizeof(N_NH));
  binFile.read(reinterpret_cast<char*>(&nHLowBound), sizeof(nHLowBound));
  binFile.read(reinterpret_cast<char*>(&nHUpBound), sizeof(nHUpBound));
  binFile.read(reinterpret_cast<char*>(&TLowBound), sizeof(TLowBound));
  binFile.read(reinterpret_cast<char*>(&TUpBound), sizeof(TUpBound));
  binFile.read(reinterpret_cast<char*>(&rLowBound), sizeof(rLowBound));
  binFile.read(reinterpret_cast<char*>(&rUpBound), sizeof(rUpBound));
  binFile.read(reinterpret_cast<char*>(&NHLowBound), sizeof(NHLowBound));
  binFile.read(reinterpret_cast<char*>(&NHUpBound), sizeof(NHUpBound));

  // Allocate memory for arrays
  float *nHGrid = new float[N_nH], *TGrid = new float[N_T], *rGrid = new float[N_r], *NHGrid = new float[N_NH];
  float *Heat = new float[N_nH * N_T * N_r * N_NH], *Cool = new float[N_nH * N_T * N_r * N_NH];
  float *muarr = new float[N_nH * N_T * N_r * N_NH], *uarr = new float[N_nH * N_T * N_r * N_NH];

  // Read grids
  binFile.read(reinterpret_cast<char*>(nHGrid), N_nH * sizeof(float));
  binFile.read(reinterpret_cast<char*>(TGrid), N_T * sizeof(float));
  binFile.read(reinterpret_cast<char*>(rGrid), N_r * sizeof(float));
  binFile.read(reinterpret_cast<char*>(NHGrid), N_NH * sizeof(float));

  // Read data arrays
  binFile.read(reinterpret_cast<char*>(Heat), N_nH * N_T * N_r * N_NH * sizeof(float));
  binFile.read(reinterpret_cast<char*>(Cool), N_nH * N_T * N_r * N_NH * sizeof(float));
  binFile.read(reinterpret_cast<char*>(muarr), N_nH * N_T * N_r * N_NH * sizeof(float));
  binFile.read(reinterpret_cast<char*>(uarr), N_nH * N_T * N_r * N_NH * sizeof(float));

  binFile.close();
  
  cout << "nHLowBound: " << nHLowBound << endl;
  cout << "nHUpBound: " << nHUpBound << endl;
  cout << "TLowBound: " << TLowBound << endl;
  cout << "TUpBound: " << TUpBound << endl;
  cout << "rLowBound: " << rLowBound << endl;
  cout << "rUpBound: " << rUpBound << endl;
  cout << "NHLowBound: " << NHLowBound << endl;
  cout << "NHUpBound: " << NHUpBound << endl << endl;
  
  //--- declaring for GPU
  float *d_nHGrid, *d_TGrid, *d_rGrid, *d_NHGrid, *d_Heat, *d_Cool, *d_muarr, *d_uarr;
  
  cudaMalloc(&d_nHGrid, N_nH * sizeof(float));
  cudaMalloc(&d_TGrid, N_T * sizeof(float));
  cudaMalloc(&d_rGrid, N_r * sizeof(float));
  cudaMalloc(&d_NHGrid, N_NH * sizeof(float));
  
  int N_tot = N_nH * N_T * N_r * N_NH;

  cudaMalloc(&d_Heat, N_tot * sizeof(float));
  cudaMalloc(&d_Cool, N_tot * sizeof(float));
  cudaMalloc(&d_muarr, N_tot * sizeof(float));
  cudaMalloc(&d_uarr, N_tot * sizeof(float));
  
  //--- Copying from Host to Device
  cudaMemcpy(d_nHGrid, nHGrid, N_nH * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_TGrid, TGrid, N_T * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rGrid, rGrid, N_r * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_NHGrid, NHGrid, N_NH * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_Heat, Heat, N_tot * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Cool, Cool, N_tot * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_muarr, muarr, N_tot * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_uarr, uarr, N_tot * sizeof(float), cudaMemcpyHostToDevice);
  //------------- End of reading and preparing HCool table file -----------
  

  //********************************************************************
  //**************** Reading the params.txt file ***********************
  //********************************************************************
  std::string filename;
  int N, ndx_BH;
  float GG, L_AGN_code_unit, M_dot_in, v_in, u_for_10K_Temp, m_sph_high_res, sigma, UnitDensity_in_cgs, Unit_u_in_cgs, unitTime_in_s,
        unitLength_in_cm;

  readParams(filename, N, ndx_BH, GG, L_AGN_code_unit, M_dot_in, v_in, u_for_10K_Temp, m_sph_high_res, sigma, UnitDensity_in_cgs, Unit_u_in_cgs, unitTime_in_s,
             unitLength_in_cm);

  std::cout << "filename: " << filename << "\n";
  std::cout << "N: " << N << "\n";
  std::cout << "ndx_BH: " << ndx_BH << "\n";
  std::cout << "GG: " << GG << "\n";
  std::cout << "L_AGN_code_unit: " << L_AGN_code_unit << "\n";
  std::cout << "M_dot_in_code_unit: " << M_dot_in << "\n";
  std::cout << "vin_in_code_unit: " << v_in << "\n";
  std::cout << "u_for_10K_Temp: " << u_for_10K_Temp << "\n";
  std::cout << "m_sph_high_res: " << m_sph_high_res << "\n";
  std::cout << "sigma: " << sigma << "\n";
  
  std::cout << "UnitDensity_in_cgs: " << UnitDensity_in_cgs << "\n";
  std::cout << "Unit_u_in_cgs: " << Unit_u_in_cgs << "\n";
  std::cout << "unitTime_in_s: " << unitTime_in_s << "\n";
  
  std::cout << "unitLength_in_cm: " << unitLength_in_cm << "\n";
  
  cout << endl;
  cout << "******************************************" << endl;
  cout << "******************************************" << endl;
  cout << "Addopted Time-step: dt_yrs = " << dt * unitTime_in_s / 3600.0 / 24.0 / 365.25 << " yrs" << endl;
  cout << "******************************************" << endl;
  cout << "******************************************" << endl;
  
  //*********************************************************************
  //******************** Reading the IC file ****************************
  //*********************************************************************
  auto data = readVectorsFromFile(filename);

  std::vector<int> &Typvec = std::get<0>(data);
  std::vector<float> &xvec = std::get<1>(data);
  std::vector<float> &yvec = std::get<2>(data);
  std::vector<float> &zvec = std::get<3>(data);
  std::vector<float> &vxvec = std::get<4>(data);
  std::vector<float> &vyvec = std::get<5>(data);
  std::vector<float> &vzvec = std::get<6>(data);
  std::vector<float> &uvec = std::get<7>(data);
  std::vector<float> &hvec = std::get<8>(data);
  std::vector<float> &epsvec = std::get<9>(data);
  std::vector<float> &massvec = std::get<10>(data);


  // declaring the arrays.
  int *Typ, *d_Typ;
  float *x, *d_x, *y, *d_y, *z, *d_z, *vx, *d_vx, *vy, *d_vy, *vz, *d_vz;
  float *mass, *d_mass, *h, *d_h, *rho, *d_rho;
  float *accx, *accy, *accz, *d_accx, *d_accy, *d_accz, *eps, *d_eps;
  float *P, *d_P, *csnd, *d_csnd, *divV, *d_divV, *curlV, *d_curlV;
  float *accx_sph, *accy_sph, *accz_sph, *d_accx_sph, *d_accy_sph, *d_accz_sph;
  float *accx_tot, *accy_tot, *accz_tot, *d_accx_tot, *d_accy_tot, *d_accz_tot;
  float *abs_acc_g, *abs_acc_tot, *v_sig, *dh_dt, *d_abs_acc_g, *d_abs_acc_tot;
  float *d_v_sig, *d_dh_dt, *u, *dudt, *d_u, *d_dudt, *utprevious;
  float *d_utprevious;
  float *Nngb_previous, *d_Nngb_previous; // Note that both are floats and not int! check smoothing func. to see why!
  float *dt_particles, *d_dt_particles;
  
  float *dudt_pre, *d_dudt_pre;
  
  float *NH, *d_NH;

  float gammah = 5.0f / 3.0f;
  float GAMMA_MINUS1 = gammah - 1.0f;


  float *gradRhoNorm, *d_gradRhoNorm;

  Typ = new int[N];

  x = new float[N];
  y = new float[N];
  z = new float[N];

  vx = new float[N];
  vy = new float[N];
  vz = new float[N];

  accx = new float[N];
  accy = new float[N];
  accz = new float[N];

  mass = new float[N];
  h = new float[N];
  rho = new float[N];
  eps = new float[N];
  P = new float[N];
  csnd = new float[N];
  
  gradRhoNorm = new float[N];

  divV = new float[N];
  curlV = new float[N];

  accx_sph = new float[N];
  accy_sph = new float[N];
  accz_sph = new float[N];

  accx_tot = new float[N];
  accy_tot = new float[N];
  accz_tot = new float[N];

  abs_acc_g = new float[N];
  abs_acc_tot = new float[N];
  v_sig = new float[N];
  dh_dt = new float[N];
  dt_particles = new float[N];

  u = new float[N];
  dudt = new float[N];
  utprevious = new float[N];
  
  dudt_pre = new float[N];

  Nngb_previous = new float[N];
  
  NH = new float[N];

  cudaMalloc(&d_Typ, N * sizeof(int));

  cudaMalloc(&d_x, N * sizeof(float));
  cudaMalloc(&d_y, N * sizeof(float));
  cudaMalloc(&d_z, N * sizeof(float));

  cudaMalloc(&d_vx, N * sizeof(float));
  cudaMalloc(&d_vy, N * sizeof(float));
  cudaMalloc(&d_vz, N * sizeof(float));

  cudaMalloc(&d_accx, N * sizeof(float));
  cudaMalloc(&d_accy, N * sizeof(float));
  cudaMalloc(&d_accz, N * sizeof(float));

  cudaMalloc(&d_mass, N * sizeof(float));
  cudaMalloc(&d_h, N * sizeof(float));
  cudaMalloc(&d_rho, N * sizeof(float));
  cudaMalloc(&d_eps, N * sizeof(float));
  cudaMalloc(&d_P, N * sizeof(float));
  cudaMalloc(&d_csnd, N * sizeof(float));
  
  cudaMalloc(&d_gradRhoNorm, N * sizeof(float));
  cudaMalloc(&d_NH, N * sizeof(float));

  cudaMalloc(&d_divV, N * sizeof(float));
  cudaMalloc(&d_curlV, N * sizeof(float));

  cudaMalloc(&d_accx_sph, N * sizeof(float));
  cudaMalloc(&d_accy_sph, N * sizeof(float));
  cudaMalloc(&d_accz_sph, N * sizeof(float));

  cudaMalloc(&d_accx_tot, N * sizeof(float));
  cudaMalloc(&d_accy_tot, N * sizeof(float));
  cudaMalloc(&d_accz_tot, N * sizeof(float));

  cudaMalloc(&d_abs_acc_g, N * sizeof(float));
  cudaMalloc(&d_abs_acc_tot, N * sizeof(float));
  cudaMalloc(&d_v_sig, N * sizeof(float));
  cudaMalloc(&d_dh_dt, N * sizeof(float));
  cudaMalloc(&d_dt_particles, N * sizeof(float));

  cudaMalloc(&d_u, N * sizeof(float));
  cudaMalloc(&d_dudt, N * sizeof(float));
  cudaMalloc(&d_utprevious, N * sizeof(float));
  
  cudaMalloc(&d_dudt_pre, N * sizeof(float));

  cudaMalloc(&d_Nngb_previous, N * sizeof(float));


  // Initialize x, y, z, etc on the Host.
  for (int i = 0; i < N; i++)
  {
    Typ[i] = Typvec[i];

    x[i] = xvec[i];
    y[i] = yvec[i];
    z[i] = zvec[i];

    vx[i] = vxvec[i];
    vy[i] = vyvec[i];
    vz[i] = vzvec[i];

    mass[i] = massvec[i];
    eps[i] = epsvec[i];

    accx[i] = 0.0f;
    accy[i] = 0.0f;
    accz[i] = 0.0f;

    accx_tot[i] = 0.0f;
    accy_tot[i] = 0.0f;
    accz_tot[i] = 0.0f;

    abs_acc_g[i] = 0.0f;
    abs_acc_tot[i] = 0.0f;
    v_sig[i] = 0.0f;

    h[i] = hvec[i]; // place holder.
    rho[i] = 0.0f;  // place holder.
    P[i] = 0.0f;    // placeholder.
    csnd[i] = 0.0f; // placeholder.
    
    gradRhoNorm[i] = 0.0f;
    NH[i] = 0.0f;

    divV[i] = 0.0f;  // placeholder.
    curlV[i] = 0.0f; // placeholder.

    accx_sph[i] = 0.0f;
    accy_sph[i] = 0.0f;
    accz_sph[i] = 0.0f;

    dh_dt[i] = 0.0f;

    u[i] = uvec[i];
    dudt[i] = 0.0f;
    utprevious[i] = 0.0f;
    
    dudt_pre[i] = 0.0f;

    dt_particles[i] = 0.0f;

    if (Typ[i] == 0)
    {
      Nngb_previous[i] = Nngb_f;
    }
    else
    {
      Nngb_previous[i] = 0.0f;
    }
  }
  

  // Copy from Host to Device.
  cudaMemcpy(d_Typ, Typ, N * sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(d_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, z, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_vx, vx, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_vy, vy, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_vz, vz, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_accx, accx, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accy, accy, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accz, accz, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_mass, mass, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_h, h, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_rho, rho, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_eps, eps, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_P, P, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_csnd, csnd, N * sizeof(float), cudaMemcpyHostToDevice);
  
  cudaMemcpy(d_gradRhoNorm, gradRhoNorm, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_NH, NH, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_divV, divV, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_curlV, curlV, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_accx_sph, accx_sph, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accy_sph, accy_sph, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accz_sph, accz_sph, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_accx_tot, accx_tot, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accy_tot, accy_tot, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accz_tot, accz_tot, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_abs_acc_g, abs_acc_g, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_abs_acc_tot, abs_acc_tot, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_v_sig, v_sig, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dh_dt, dh_dt, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dt_particles, dt_particles, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_u, u, N * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dudt, dudt, N * sizeof(float), cudaMemcpyHostToDevice);
  
  cudaMemcpy(d_utprevious, utprevious, N * sizeof(float), cudaMemcpyHostToDevice);
  
  cudaMemcpy(d_dudt_pre, dudt_pre, N * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_Nngb_previous, Nngb_previous, N * sizeof(float), cudaMemcpyHostToDevice);
  
  
  // allocate device data
  float *d_leftX;
  float *d_rightX;
  float *d_bottomX;
  float *d_topX;
  float *d_frontX;
  float *d_backX;
  cudaMalloc((void**)&d_leftX, sizeof(float));
  cudaMalloc((void**)&d_rightX, sizeof(float));
  cudaMalloc((void**)&d_bottomX, sizeof(float));
  cudaMalloc((void**)&d_topX, sizeof(float));
  cudaMalloc((void**)&d_frontX, sizeof(float));
  cudaMalloc((void**)&d_backX, sizeof(float));

  cudaMemset(d_leftX, 0, sizeof(float));
  cudaMemset(d_rightX, 0, sizeof(float));
  cudaMemset(d_bottomX, 0, sizeof(float));
  cudaMemset(d_topX, 0, sizeof(float));
  cudaMemset(d_frontX, 0, sizeof(float));
  cudaMemset(d_backX, 0, sizeof(float));
  
  
  
  //============== For BH Tree =================
  int numParticles, numNodes;

  float *h_mass;
  float *h_x;
  float *h_y;
  float *h_z;
  float *h_ax;
  float *h_ay;
  float *h_az;

  int *h_child;
  int *h_start;
  int *h_sorted;
  int *h_count;

  float *d_left;
  float *d_right;
  float *d_bottom;
  float *d_top;
  float *d_front;
  float *d_back;

  float *dev_mass;
  float *dev_x;
  float *dev_y;
  float *dev_z;
  float *dev_ax;
  float *dev_ay;
  float *dev_az;

  int *d_index;
  int *d_child;
  int *d_start;
  int *d_sorted;
  int *d_count;

  int *d_mutex;  //used for locking
  
  int blockSize_bh, gridSize_bh;
  int nBodies;
  //===================
    
  

  //int blockSize = 256;                            // number of threads in a block
  int gridSize = (N + blockSize - 1) / blockSize; // Number of blocks in a grid

  const float visc_alpha = 2.0f;

  float t;

  t = 0.0f;

  float tEnd = 1.0f;
  float Nt = ceil(tEnd / dt) + 1;

  //float Zmetal = 0.1; // ==> [Z/H] = -1.
  
  
  //----- Used for the possible need to decrease the value of nSplit -------
  int flag;
  int *dev_flag;
  cudaMalloc((void **)&dev_flag, sizeof(int));
  cudaMemset(dev_flag, 0, sizeof(int)); // Initialize to 0
  //-------------------------------------------------------------------------
  
  
  //------- For ngb_new -----
  //--- Finding max_h
  float max_h = h[0];
  for (int i = 1; i < N; i++)
  {
    if (h[i] > max_h)
      max_h = h[i];
  }
  
  float binSize = 7.0f; // For Lim = 1 ==> binSize = 3. for Lim = 2 ==> binSize = 5. For Lim = 3 ==> binSize = 7
  float hcoeff = 2.0f; // The higher the better. But 3 should be fine. Higher values decrease nSplit leading to slower code!
  //----- Finding L_box
  float max_x = 0.0;
  float max_y = 0.0;
  float max_z = 0.0;
  for (int i = 1; i < N; i++)
  {
    if (abs(x[i]) > max_x)
      max_x = abs(x[i]);
    if (abs(y[i]) > max_y)
      max_y = abs(y[i]);
    if (abs(z[i]) > max_z)
      max_z = abs(z[i]);
  }
  float L_box = max(max(max_x, max_y), max_z);
  
  int nSplit = binSize * L_box / pow(hcoeff, 1.0f/3.0f) / 2.0f / max_h;
  cout << "Adopted nSplit = " << nSplit << endl;
  
  //int nSplit = 40; //32; // for 3D, never go above 10 as 10^3 = 1000 and the maximum blockSize in GPUs is 1024!! IMPORTANT !!
  int Ncell = nSplit * nSplit * nSplit;
  
  //---- Finding nBodies (i.e. Those with Typ == 0 EXCLUDING BHole) REPEAT in each time-step -----
  int nBodiesB = 0;
  for (int i = 0; i < N; i++)
  {
    if (Typ[i] != -1) // !!!!!!!!!!BH should be excluded !!!===> No we can include BH in this step as will be excluded in smoothing_ step automatically!!
      nBodiesB++;
  }
  
  int *groupedIndex, *d_groupedIndex, *countx, *d_countx;
  
  groupedIndex = new int[nBodiesB];
  countx = new int[2*Ncell];
  
  cudaMalloc((void **)&d_groupedIndex, nBodiesB * sizeof(int));
  cudaMalloc((void **)&d_countx, 2*Ncell * sizeof(int));

  for (int i = 0; i < nBodiesB; i++)
  {
    groupedIndex[i] = 0; // just a place holder. Its initial values don't matter as it will be replaced anyway!
  }
  
  for (int i = 0; i < (2*Ncell); i++)
  {
    countx[i] = 0;
  }
  
  cudaMemcpy(d_groupedIndex, groupedIndex, nBodiesB * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_countx, countx, 2*Ncell * sizeof(int), cudaMemcpyHostToDevice);


  //----- computing maxCoord ------
  cudaMemcpy(Typ, d_Typ, N * sizeof(int), cudaMemcpyDeviceToHost);
  nBodies = 0;
  for (int i = 0; i < N; i++)
  {
    if (Typ[i] != -1)
    {
      nBodies++;
    }
  }
  int *d_mutexX;
  cudaMalloc((void**)&d_mutexX, sizeof(int));
  gpuErrchk(cudaMemset(d_mutexX, 0, sizeof(int)));
  compute_bounding_box_kernel<<< 1, blockSize >>>(d_mutexX, d_x, d_y, d_z, d_leftX, d_rightX, d_bottomX, d_topX, d_frontX, d_backX,
                                                  nBodies);
  cudaDeviceSynchronize();
  
  float *h_leftX;
  float *h_rightX;
  float *h_bottomX;
  float *h_topX;
  float *h_frontX;
  float *h_backX;
  // allocate host data
  h_leftX = new float;
  h_rightX = new float;
  h_bottomX = new float;
  h_topX = new float;
  h_frontX = new float;
  h_backX = new float;
  
  cudaMemcpy(h_leftX, d_leftX, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_rightX, d_rightX, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_bottomX, d_bottomX, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_topX, d_topX, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_frontX, d_frontX, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_backX, d_backX, sizeof(float), cudaMemcpyDeviceToHost);
  
  printf("\n");
  printf("h_left, h_right, h_bottom, h_top, h_front, h_back = %f, %f, %f, %f, %f, %f\n",h_leftX[0],h_rightX[0], h_bottomX[0], h_topX[0], h_frontX[0], h_backX[0]);
  printf("\n");
  printf("\n");
  
  float maxRange = fmax(abs(*h_leftX), fmax(abs(*h_rightX), fmax(abs(*h_bottomX), fmax(abs(*h_topX), fmax(abs(*h_frontX), abs(*h_backX))))));
  float W_cell = ceil(2.0 * maxRange) / nSplit;
  
  float x_min = -1.0 * maxRange;
  float y_min = -1.0 * maxRange;
  float z_min = -1.0 * maxRange;
  
  cout << "x_min, y_min, z_min = " << x_min << endl;
  
  CountBodies<<<1, 512>>>(d_x, d_y, d_z, x_min, y_min, z_min, W_cell, nSplit, d_countx, Ncell, nBodiesB);
  cudaDeviceSynchronize();
  
    
  int blockSize_CO = 512;
  int gridSize_CO = (Ncell + blockSize_CO - 1) / blockSize_CO;
  ComputeOffset<<<gridSize_CO, blockSize_CO>>>(d_countx, Ncell);
  cudaDeviceSynchronize();
  
  GroupBodies<<<1, 1024>>>(d_x, d_y, d_z, d_groupedIndex, x_min, y_min, z_min, W_cell, nSplit, d_countx, Ncell, nBodiesB);
  cudaDeviceSynchronize();

  
  cudaMemcpy(countx, d_countx, 2*Ncell * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(groupedIndex, d_groupedIndex, nBodiesB * sizeof(int), cudaMemcpyDeviceToHost);
  
  int *offSet;
  offSet = new int[Ncell+1];
  
  offSet[0] = 0;
  for (int i = Ncell; i < 2*Ncell; i++)
  {
    offSet[i-Ncell+1] = countx[i];
  }

  int *d_offSet;
  cudaMalloc((void **)&d_offSet, (Ncell+1) * sizeof(int));
  cudaMemcpy(d_offSet, offSet, (Ncell+1) * sizeof(int), cudaMemcpyHostToDevice);
  //---------------------------

  //----------
  int MAX_ngb = 200;
  int MAX_N = N * MAX_ngb;
  int *ngb, *d_ngb;
  
  ngb = new int[MAX_N];
  
  cudaMalloc(&d_ngb, MAX_N * sizeof(int));
  
  for (int i = 0; i < MAX_N; i++)
  {
    ngb[i] = -1;
  }

  cudaMemcpy(d_ngb, ngb, MAX_N * sizeof(int), cudaMemcpyHostToDevice);
  //----------

  //------------------------------------------------------
  //------------------- ngbDB_new_v3 ---------------------
  //------------------------------------------------------
  auto T_ngb_new = std::chrono::high_resolution_clock::now();
  ngbDB_new_v3<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_h,
                                        x_min, y_min, z_min, W_cell, nSplit,
                                        d_offSet, d_groupedIndex,
                                        d_ngb, MAX_ngb, N, dev_flag);
  cudaDeviceSynchronize();
  auto end_ngb_new = std::chrono::high_resolution_clock::now();
  auto elapsed_ngb_new = std::chrono::duration_cast<std::chrono::nanoseconds>(end_ngb_new - T_ngb_new);
  cout << "T_ngb_new = " << elapsed_ngb_new.count() * 1e-9 << endl;

  cudaMemcpy(&flag, dev_flag, sizeof(int), cudaMemcpyDeviceToHost);
  if (flag == 1)
  {
    std::cerr << "Exiting program: nSplit value is large - decrease it !!!" << std::endl;
    cudaFree(dev_flag);
    exit(EXIT_FAILURE);
  }

  //-----------------------------------------------
  //-------------- Smoothing Length_ngb ---------------
  //-----------------------------------------------
  auto T_hh = std::chrono::high_resolution_clock::now();
  smoothing_h_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_h,
                                           N, Ndown, Nup, coeff,
                                           Nngb_f, d_Nngb_previous, d_divV, d_ngb, MAX_ngb, dt);
  cudaDeviceSynchronize();
  auto end_hh = std::chrono::high_resolution_clock::now();
  auto elapsed_hh = std::chrono::duration_cast<std::chrono::nanoseconds>(end_hh - T_hh);
  cout << "T_h = " << elapsed_hh.count() * 1e-9 << endl;
  
  //-----------------------------------------------
  //----------------- getDensity ------------------
  //-----------------------------------------------
  getDensity_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_mass,
                                          d_rho, d_h, d_ngb, MAX_ngb, N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //---------------- getPressure ------------------
  //-----------------------------------------------
  getPressure_Adiabatic<<<gridSize, blockSize>>>(d_Typ, d_P, d_rho, d_u, gammah, N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //----------------- getCsound -------------------
  //-----------------------------------------------
  getCsound_Adiabatic<<<gridSize, blockSize>>>(d_Typ, d_csnd, d_u, gammah, N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //----------------- div_curlV -------------------
  //-----------------------------------------------
  div_curlVel_ngb<<<gridSize, blockSize>>>(d_Typ, d_divV, d_curlV, d_x, d_y, d_z, d_vx, d_vy, d_vz,
                                           d_rho, d_mass, d_ngb, MAX_ngb, d_h, N);
  cudaDeviceSynchronize();
  
  //-----------------------------------------------
  //------------------ acc_sph --------------------
  //-----------------------------------------------
  acc_sph_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                       d_divV, d_curlV, d_mass, d_P, d_accx_sph, d_accy_sph,
                                       d_accz_sph, d_ngb, MAX_ngb, visc_alpha, N);
  cudaDeviceSynchronize();


  //-----------------------------------------------
  //------------------ getAcc_g -------------------
  //-----------------------------------------------
  cudaMemcpy(Typ, d_Typ, N * sizeof(int), cudaMemcpyDeviceToHost);
  nBodies = 0;
  for (int i = 0; i < N; i++)
  {
    if (Typ[i] != -1)
    {
      nBodies++;
    }
  }
  
  cout << "nBodies = " << nBodies << endl;
  
  numParticles = nBodies; // nBodies is the number of patticles with Typ != -1.
  int eXtraSpace = 10000000;
  numNodes = 2 * numParticles + eXtraSpace;

  blockSize_bh = blockSize;
  gridSize_bh = (numParticles + blockSize_bh - 1) / blockSize_bh;
  
  // allocate host data
  h_mass = new float[numNodes];
  h_x = new float[numNodes];
  h_y = new float[numNodes];
  h_z = new float[numNodes];
  h_ax = new float[numNodes];
  h_ay = new float[numNodes];
  h_az = new float[numNodes];
  h_child = new int[8*numNodes];
  h_start = new int[numNodes];
  h_sorted = new int[numNodes];
  h_count = new int[numNodes];

  // allocate device data
  gpuErrchk(cudaMalloc((void**)&d_left, sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&d_right, sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&d_bottom, sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&d_top, sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&d_front, sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&d_back, sizeof(float)));

  gpuErrchk(cudaMemset(d_left, 0, sizeof(float)));
  gpuErrchk(cudaMemset(d_right, 0, sizeof(float)));
  gpuErrchk(cudaMemset(d_bottom, 0, sizeof(float)));
  gpuErrchk(cudaMemset(d_top, 0, sizeof(float)));
  gpuErrchk(cudaMemset(d_front, 0, sizeof(float)));
  gpuErrchk(cudaMemset(d_back, 0, sizeof(float)));

  gpuErrchk(cudaMalloc((void**)&dev_mass, numNodes*sizeof(float)));

  gpuErrchk(cudaMalloc((void**)&dev_x, numNodes*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&dev_y, numNodes*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&dev_z, numNodes*sizeof(float)));

  gpuErrchk(cudaMalloc((void**)&dev_ax, numNodes*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&dev_ay, numNodes*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&dev_az, numNodes*sizeof(float)));

  gpuErrchk(cudaMalloc((void**)&d_index, sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&d_child, 8*numNodes*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&d_start, numNodes*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&d_sorted, numNodes*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&d_count, numNodes*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&d_mutex, sizeof(int))); 

  gpuErrchk(cudaMemset(d_start, -1, numNodes*sizeof(int)));
  gpuErrchk(cudaMemset(d_sorted, 0, numNodes*sizeof(int)));
  gpuErrchk(cudaMemset(d_mutex, 0, sizeof(int)));
 
  reset_arrays_kernel<<< 1, blockSize_bh >>>(d_mutex, dev_x, dev_y, dev_z, dev_mass, d_count, d_start, d_sorted, d_child, d_index,
                                             d_left, d_right, d_bottom, d_top, d_front, d_back, numParticles, numNodes);
  cudaDeviceSynchronize();
  
  // initializing x, y, z, mass -----
  cudaMemcpy(x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(y, d_y, N * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(z, d_z, N * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(mass, d_mass, N * sizeof(float), cudaMemcpyDeviceToHost);
    
  for (int i = 0; i < numParticles; i++)
  {
    h_x[i] = x[i];
    h_y[i] = y[i];
    h_z[i] = z[i];
    
    h_mass[i] = mass[i];
  }

  cudaMemcpy(dev_x, h_x, numNodes * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y, h_y, numNodes * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z, h_z, numNodes * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_mass, h_mass, numNodes * sizeof(float), cudaMemcpyHostToDevice);
  
  compute_bounding_box_kernel<<< 1, blockSize_bh >>>(d_mutex, dev_x, dev_y, dev_z, d_left, d_right, d_bottom, d_top, d_front, d_back,
                                                     numParticles);
  cudaDeviceSynchronize();
  
  float *h_left;
  float *h_right;
  float *h_bottom;
  float *h_top;
  float *h_front;
  float *h_back;
  // allocate host data
  h_left = new float;
  h_right = new float;
  h_bottom = new float;
  h_top = new float;
  h_front = new float;
  h_back = new float;
  
  cudaMemcpy(h_left, d_left, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_right, d_right, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_bottom, d_bottom, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_top, d_top, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_front, d_front, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_back, d_back, sizeof(float), cudaMemcpyDeviceToHost);
  
  int *h_index = new int;
  cudaMemcpy(h_index, d_index, sizeof(int), cudaMemcpyDeviceToHost);
  printf("\n");
  printf("h_left, h_right, h_bottom, h_top, h_front, h_back = %f, %f, %f, %f, %f, %f\n", h_left[0], h_right[0], h_bottom[0], h_top[0], h_front[0], h_back[0]);
  printf("\n");
  printf("initial index = %d\n", h_index[0]);
  printf("\n");

  auto T_build_tree_kernel = std::chrono::high_resolution_clock::now();
  build_tree_kernel<<< 1, 256 >>>(dev_x, dev_y, dev_z, dev_mass, d_count, d_start, d_child, d_index, d_left, d_right, d_bottom, d_top, d_front, d_back,
                                  numParticles, numNodes);
  cudaDeviceSynchronize();  
  auto end_build_tree_kernel = std::chrono::high_resolution_clock::now();
  auto elapsed_build_tree_kernel = std::chrono::duration_cast<std::chrono::nanoseconds>(end_build_tree_kernel - T_build_tree_kernel);
  cout << "Tree construction time = " << elapsed_build_tree_kernel.count() * 1e-9 << endl;

  centre_of_mass_kernel<<<1, blockSize_bh>>>(dev_x, dev_y, dev_z, dev_mass, d_index, numParticles);
  cudaDeviceSynchronize();  
  
  sort_kernel<<< 1, 256 >>>(d_count, d_start, d_sorted, d_child, d_index, numParticles);
  cudaDeviceSynchronize();
  
  auto T_Force = std::chrono::high_resolution_clock::now();
  compute_forces_kernel<<< gridSize_bh, blockSize_bh >>>(dev_x, dev_y, dev_z, dev_ax, dev_ay, dev_az, dev_mass, d_eps, d_sorted, d_child,
                                                         d_left, d_right, d_bottom, d_top, d_front, d_back, numParticles);
  cudaDeviceSynchronize();
  auto end_Force = std::chrono::high_resolution_clock::now();
  auto elapsed_Force = std::chrono::duration_cast<std::chrono::nanoseconds>(end_Force - T_Force);
  cout << "T_Force = " << elapsed_Force.count() * 1e-9 << endl;
  
  cudaMemcpy(h_ax, dev_ax, numNodes * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_ay, dev_ay, numNodes * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_az, dev_az, numNodes * sizeof(float), cudaMemcpyDeviceToHost);

  for (int i = 0; i < numParticles; i++)
  {
    accx[i] = h_ax[i];
    accy[i] = h_ay[i];
    accz[i] = h_az[i];
  }
  
  cudaMemcpy(d_accx, accx, numParticles * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accy, accy, numParticles * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_accz, accz, numParticles * sizeof(float), cudaMemcpyHostToDevice);
  
  delete[] h_mass;
  delete[] h_x;
  delete[] h_y;
  delete[] h_z;
  delete[] h_ax;
  delete[] h_ay;
  delete[] h_az;
  delete[] h_child;
  delete[] h_start;
  delete[] h_sorted;
  delete[] h_count;
  
  cudaFree(d_left);
  cudaFree(d_right);
  cudaFree(d_bottom);
  cudaFree(d_top);
  cudaFree(d_front);
  cudaFree(d_back);
  
  cudaFree(dev_mass);
  cudaFree(dev_x);
  cudaFree(dev_y);
  cudaFree(dev_z);
  
  cudaFree(dev_ax);
  cudaFree(dev_ay);
  cudaFree(dev_az);
  
  cudaFree(d_index);
  cudaFree(d_child);
  cudaFree(d_start);
  cudaFree(d_sorted);
  cudaFree(d_count);
  cudaFree(d_mutex);

  //-----------------------------------------------
  //------------------ acc_tot --------------------
  //-----------------------------------------------
  acc_g_sph<<<gridSize, blockSize>>>(d_Typ, d_accx_tot, d_accy_tot, d_accz_tot,
                                     d_accx, d_accy, d_accz,
                                     d_accx_sph, d_accy_sph, d_accz_sph,
                                     N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //------------------- du_dt ---------------------
  //-----------------------------------------------  
  get_dU_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                      d_divV, d_curlV, d_mass, d_P, d_dudt,
                                      d_ngb, MAX_ngb, visc_alpha, N);
  cudaDeviceSynchronize();
  
  

  //-----------------------------------------------
  //---------------- u evolution ------------------
  //-----------------------------------------------

  u_updater<<<gridSize, blockSize>>>(d_Typ, d_u, d_dudt, d_utprevious, dt, N);
  cudaDeviceSynchronize();

  //const float C_CFL = 0.25;
  
  delete[] ngb;
  cudaFree(d_ngb);
  
  delete[] groupedIndex;
  delete[] countx;
  delete[] offSet;
  
  cudaFree(d_groupedIndex);
  cudaFree(d_countx);
  cudaFree(d_offSet);

  float leftover_mass = 0.0f;
  float *d_leftover_mass;
  cudaMalloc((void **)&d_leftover_mass, sizeof(float));
  cudaMemcpy(d_leftover_mass, &leftover_mass, sizeof(float), cudaMemcpyHostToDevice);

  //float RR0 = maxRange - 0.02; // particles at this radius will not be moved.... my boundary condition to avoid very low density and nsplit problem!!!

  // **************************************************************
  // *********************** MAIN LOOP ****************************
  // **************************************************************

  int counter = 0; // This is used to save fewer output files, e.g. 1 snap-shot per 2 time-step!

  while (t < tEnd)
  {
  
    auto begin = std::chrono::high_resolution_clock::now();
  
    //----------
    int *ngb, *d_ngb;
    ngb = new int[MAX_N];

    cudaMalloc(&d_ngb, MAX_N * sizeof(int));

    auto T_MovingData_1 = std::chrono::high_resolution_clock::now();
    cudaMemset(d_ngb, -1, MAX_N * sizeof(int));
    auto end_MovingData_1 = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData_1 - T_MovingData_1);
    cout << "T_MovingData_1 = " << elapsed_MovingData_1.count() * 1e-9 << endl;
    //----------

    //****************** velocity evolution *******************
    v_evolve<<<gridSize, blockSize>>>(d_Typ, d_vx, d_vy, d_vz, d_accx_tot, d_accy_tot,
                                      d_accz_tot, dt, N);
    cudaDeviceSynchronize();

    //****************** position evolution (BH fixed at [0, 0, 0]) *******************
    
    r_evolve<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, dt, ndx_BH, N);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //r_evolve_fixPos_atRadiusRR0<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, dt, ndx_BH, N, RR0); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cudaDeviceSynchronize(); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    //------- For ngb_new -----
    auto T_nSplit = std::chrono::high_resolution_clock::now();
    //--- Finding max_h
    max_h = h[0];
    for (int i = 1; i < N; i++)
    {
      if (h[i] > max_h)
        max_h = h[i];
    }
    
    L_box = maxRange;
    
    int nSplit = binSize * L_box / pow(hcoeff, 1.0f/3.0f) / 2.0f / max_h;
    cout << "Adopted nSplit = " << nSplit << endl;
    Ncell = nSplit * nSplit * nSplit;
    auto end_nSplit = std::chrono::high_resolution_clock::now();
    auto elapsed_nSplit = std::chrono::duration_cast<std::chrono::nanoseconds>(end_nSplit - T_nSplit);
    cout << "T_nSplit = " << elapsed_nSplit.count() * 1e-9 << endl;
    
    //---- Finding nBodies (i.e. Those with Typ == 0 EXCLUDING BHole) REPEAT in each time-step -----
    nBodiesB = 0;
    for (int i = 0; i < N; i++)
    {
      if (Typ[i] != -1) // BH should be excluded !!!
        nBodiesB++;
    }
    
    groupedIndex = new int[nBodiesB];
    countx = new int[2*Ncell];
    
    cudaMalloc((void **)&d_groupedIndex, nBodiesB * sizeof(int));
    cudaMalloc((void **)&d_countx, 2*Ncell * sizeof(int));

    for (int i = 0; i < nBodiesB; i++)
    {
      groupedIndex[i] = 0; // just a place holder. Its initial values don't matter as it will be replaced anyway!
    }
    
    for (int i = 0; i < (2*Ncell); i++)
    {
      countx[i] = 0;
    }
    
    auto T_MovingData_3 = std::chrono::high_resolution_clock::now();
    cudaMemcpy(d_groupedIndex, groupedIndex, nBodiesB * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_countx, countx, 2*Ncell * sizeof(int), cudaMemcpyHostToDevice);
    auto end_MovingData_3 = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData_3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData_3 - T_MovingData_3);
    cout << "T_MovingData_3 = " << elapsed_MovingData_3.count() * 1e-9 << endl;

    //----- computing maxCoord ------
    cudaMemset(d_leftX, 0, sizeof(float));
    cudaMemset(d_rightX, 0, sizeof(float));
    cudaMemset(d_bottomX, 0, sizeof(float));
    cudaMemset(d_topX, 0, sizeof(float));
    cudaMemset(d_frontX, 0, sizeof(float));
    cudaMemset(d_backX, 0, sizeof(float));
    
    cudaMemcpy(Typ, d_Typ, N * sizeof(int), cudaMemcpyDeviceToHost);
    nBodies = 0;
    for (int i = 0; i < N; i++)
    {
      if (Typ[i] != -1)
      {
        nBodies++;
      }
    }

    gpuErrchk(cudaMemset(d_mutexX, 0, sizeof(int)));
    compute_bounding_box_kernel<<< 1, blockSize >>>(d_mutexX, d_x, d_y, d_z, d_leftX, d_rightX, d_bottomX, d_topX, d_frontX, d_backX,
                                                    nBodies);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_leftX, d_leftX, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_rightX, d_rightX, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_bottomX, d_bottomX, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_topX, d_topX, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_frontX, d_frontX, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_backX, d_backX, sizeof(float), cudaMemcpyDeviceToHost);
    
    float maxRange = fmax(abs(*h_leftX), fmax(abs(*h_rightX), fmax(abs(*h_bottomX), fmax(abs(*h_topX), fmax(abs(*h_frontX), abs(*h_backX))))));
    W_cell = ceil(2.0 * maxRange) / nSplit;
    
    x_min = -1.0 * maxRange;
    y_min = -1.0 * maxRange;
    z_min = -1.0 * maxRange;
    
    CountBodies<<<1, 512>>>(d_x, d_y, d_z, x_min, y_min, z_min, W_cell, nSplit, d_countx, Ncell, nBodiesB);
    cudaDeviceSynchronize();
    
      
    blockSize_CO = 512;
    gridSize_CO = (Ncell + blockSize_CO - 1) / blockSize_CO;
    ComputeOffset<<<gridSize_CO, blockSize_CO>>>(d_countx, Ncell);
    cudaDeviceSynchronize();
    
    GroupBodies<<<1, 1024>>>(d_x, d_y, d_z, d_groupedIndex, x_min, y_min, z_min, W_cell, nSplit, d_countx, Ncell, nBodiesB);
    cudaDeviceSynchronize();

    
    cudaMemcpy(countx, d_countx, 2*Ncell * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(groupedIndex, d_groupedIndex, nBodiesB * sizeof(int), cudaMemcpyDeviceToHost);
    
    int *offSet;
    offSet = new int[Ncell+1];
    
    offSet[0] = 0;
    for (int i = Ncell; i < 2*Ncell; i++)
    {
      offSet[i-Ncell+1] = countx[i];
    }
    
    int *d_offSet;
    cudaMalloc((void **)&d_offSet, (Ncell+1) * sizeof(int));
    cudaMemcpy(d_offSet, offSet, (Ncell+1) * sizeof(int), cudaMemcpyHostToDevice);
    //---------------------------

    //------------------------------------------------------
    //------------------- ngbDB_new_v3 ---------------------
    //------------------------------------------------------
    auto T_ngb_new = std::chrono::high_resolution_clock::now();
    ngbDB_new_v3<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_h,
                                          x_min, y_min, z_min, W_cell, nSplit,
                                          d_offSet, d_groupedIndex,
                                          d_ngb, MAX_ngb, N, dev_flag);
    cudaDeviceSynchronize();
    auto end_ngb_new = std::chrono::high_resolution_clock::now();
    auto elapsed_ngb_new = std::chrono::duration_cast<std::chrono::nanoseconds>(end_ngb_new - T_ngb_new);
    cout << "T_ngb_new = " << elapsed_ngb_new.count() * 1e-9 << endl;

    cudaMemcpy(&flag, dev_flag, sizeof(int), cudaMemcpyDeviceToHost);
    /*
    if (flag == 1)
    {
      std::cerr << "Exiting program: nSplit value is large - decrease it !!!" << std::endl;
      cudaFree(dev_flag);
      exit(EXIT_FAILURE);
    }
    */
 

    //****************** Smoothing Length *********************
    auto T_hh = std::chrono::high_resolution_clock::now();
    smoothing_h_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_h,
                                             N, Ndown, Nup, coeff,
                                             Nngb_f, d_Nngb_previous, d_divV, d_ngb, MAX_ngb, dt);
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) 
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cout << "Smoothing Length ----> Successful!." << endl;
    auto end_hh = std::chrono::high_resolution_clock::now();
    auto elapsed_hh = std::chrono::duration_cast<std::chrono::nanoseconds>(end_hh - T_hh);
    cout << "T_h = " << elapsed_hh.count() * 1e-9 << endl;




    //****************** Set eps of Gas equal to h ******************
    set_eps_of_gas_to_h<<<gridSize, blockSize>>>(d_Typ, d_eps, d_h, N);
    cudaDeviceSynchronize();

    //****************** getDensity ***********************
    auto T_density = std::chrono::high_resolution_clock::now();
    getDensity_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_mass,
                                            d_rho, d_h, d_ngb, MAX_ngb, N);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) 
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cout << "getDensity ----> Successful!." << endl;
    auto end_density = std::chrono::high_resolution_clock::now();
    auto elapsed_density = std::chrono::duration_cast<std::chrono::nanoseconds>(end_density - T_density);
    cout << "T_density = " << elapsed_density.count() * 1e-9 << endl;
    
    
    
    
    //****************** grad_rho_norm_ngb ***********************
    auto T_grad_rho = std::chrono::high_resolution_clock::now();
    grad_rho_norm_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_rho, d_mass,
                                               d_ngb, MAX_ngb, d_h, d_gradRhoNorm, N);
    cudaDeviceSynchronize();
    
    cudaMemcpy(rho, d_rho, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(gradRhoNorm, d_gradRhoNorm, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h, d_h, N * sizeof(float), cudaMemcpyDeviceToHost);
    
    float NH_tmp = 0.0, h_inter = 0.0, nH_tmp = 0.0;
    
    for (int i = 0; i < N; i++)
    {
      nH_tmp = rho[i] * UnitDensity_in_cgs * XH / mH;
    
      h_inter = h[i] / cbrtf(Nngb_f); // cbrtf is CUDA function to get cube root for single-precision floating-point numbers!
      
      NH_tmp = 0.5 * (h_inter + rho[i] / gradRhoNorm[i]) * unitLength_in_cm * nH_tmp;
      NH[i] = log10(NH_tmp);
    }
    
    cudaMemcpy(d_NH, NH, N * sizeof(float), cudaMemcpyHostToDevice);
    auto end_grad_rho = std::chrono::high_resolution_clock::now();
    auto elapsed_grad_rho = std::chrono::duration_cast<std::chrono::nanoseconds>(end_grad_rho - T_grad_rho);
    cout << "T_grad_rho = " << elapsed_grad_rho.count() * 1e-9 << endl;

    cout << endl;
    for (int i = 0; i < 4; i++)
    {
      cout << "NH = " << NH[i] << endl;
    }
    cout << endl;
    
    
    
    
    
    
    //****************** getPressure **********************
    getPressure_Adiabatic<<<gridSize, blockSize>>>(d_Typ, d_P, d_rho, d_u, gammah, N);
    cudaDeviceSynchronize();

    //****************** getCsound ************************
    getCsound_Adiabatic<<<gridSize, blockSize>>>(d_Typ, d_csnd, d_u, gammah, N);
    cudaDeviceSynchronize();

    //****************** div_curlVX ************************
    auto T_divCurl = std::chrono::high_resolution_clock::now();
    div_curlVel_ngb<<<gridSize, blockSize>>>(d_Typ, d_divV, d_curlV, d_x, d_y, d_z, d_vx, d_vy, d_vz,
                                             d_rho, d_mass, d_ngb, MAX_ngb, d_h, N);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) 
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cout << "divV ----> Successful!." << endl;
    auto end_divCurl = std::chrono::high_resolution_clock::now();
    auto elapsed_divCurl = std::chrono::duration_cast<std::chrono::nanoseconds>(end_divCurl - T_divCurl);
    cout << "T_divCurl = " << elapsed_divCurl.count() * 1e-9 << endl;
    
    //****************** acc_sphX **************************
    auto T_acc_sph = std::chrono::high_resolution_clock::now();
    acc_sph_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                         d_divV, d_curlV, d_mass, d_P, d_accx_sph, d_accy_sph,
                                         d_accz_sph, d_ngb, MAX_ngb, visc_alpha, N);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) 
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cout << "acc_sph ----> Successful!." << endl;
    auto end_acc_sph = std::chrono::high_resolution_clock::now();
    auto elapsed_acc_sph = std::chrono::duration_cast<std::chrono::nanoseconds>(end_acc_sph - T_acc_sph);
    cout << "T_acc_sph = " << elapsed_acc_sph.count() * 1e-9 << endl;
    
    
    
    //****************** getAcc_gX *************************
    auto T_acc_g_Full_process = std::chrono::high_resolution_clock::now();
    cudaMemcpy(Typ, d_Typ, N * sizeof(int), cudaMemcpyDeviceToHost);
    nBodies = 0;
    for (int i = 0; i < N; i++)
    {
      if (Typ[i] != -1)
      {
        nBodies++;
      }
    }
    
    cout << "nBodies = " << nBodies << endl;
    
    numParticles = nBodies; // nBodies is the number of patticles with Typ != -1.
    numNodes = 2 * numParticles + eXtraSpace;

    blockSize_bh = blockSize;
    gridSize_bh = (numParticles + blockSize_bh - 1) / blockSize_bh;
    
    // allocate host data
    h_mass = new float[numNodes];
    h_x = new float[numNodes];
    h_y = new float[numNodes];
    h_z = new float[numNodes];
    h_ax = new float[numNodes];
    h_ay = new float[numNodes];
    h_az = new float[numNodes];
    h_child = new int[8*numNodes];
    h_start = new int[numNodes];
    h_sorted = new int[numNodes];
    h_count = new int[numNodes];

    // allocate device data
    gpuErrchk(cudaMalloc((void**)&d_left, sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&d_right, sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&d_bottom, sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&d_top, sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&d_front, sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&d_back, sizeof(float)));

    gpuErrchk(cudaMemset(d_left, 0, sizeof(float)));
    gpuErrchk(cudaMemset(d_right, 0, sizeof(float)));
    gpuErrchk(cudaMemset(d_bottom, 0, sizeof(float)));
    gpuErrchk(cudaMemset(d_top, 0, sizeof(float)));
    gpuErrchk(cudaMemset(d_front, 0, sizeof(float)));
    gpuErrchk(cudaMemset(d_back, 0, sizeof(float)));

    gpuErrchk(cudaMalloc((void**)&dev_mass, numNodes*sizeof(float)));

    gpuErrchk(cudaMalloc((void**)&dev_x, numNodes*sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&dev_y, numNodes*sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&dev_z, numNodes*sizeof(float)));

    gpuErrchk(cudaMalloc((void**)&dev_ax, numNodes*sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&dev_ay, numNodes*sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&dev_az, numNodes*sizeof(float)));

    gpuErrchk(cudaMalloc((void**)&d_index, sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&d_child, 8*numNodes*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&d_start, numNodes*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&d_sorted, numNodes*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&d_count, numNodes*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&d_mutex, sizeof(int))); 

    gpuErrchk(cudaMemset(d_start, -1, numNodes*sizeof(int)));
    gpuErrchk(cudaMemset(d_sorted, 0, numNodes*sizeof(int)));
    gpuErrchk(cudaMemset(d_mutex, 0, sizeof(int)));
   
    auto T_reset_kernel = std::chrono::high_resolution_clock::now();
    reset_arrays_kernel<<< 1, blockSize_bh >>>(d_mutex, dev_x, dev_y, dev_z, dev_mass, d_count, d_start, d_sorted, d_child, d_index,
                                               d_left, d_right, d_bottom, d_top, d_front, d_back, numParticles, numNodes);
    cudaDeviceSynchronize();
    auto end_reset_kernel = std::chrono::high_resolution_clock::now();
    auto elapsed_reset_kernel = std::chrono::duration_cast<std::chrono::nanoseconds>(end_reset_kernel - T_reset_kernel);
    cout << "T_reset_kernel = " << elapsed_reset_kernel.count() * 1e-9 << endl;
    
    // initializing x, y, z, mass -----
    auto T_MovingData_5 = std::chrono::high_resolution_clock::now();
    cudaMemcpy(x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(y, d_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(z, d_z, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(mass, d_mass, N * sizeof(float), cudaMemcpyDeviceToHost);
    auto end_MovingData_5 = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData_5 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData_5 - T_MovingData_5);
    cout << "T_MovingData_5 = " << elapsed_MovingData_5.count() * 1e-9 << endl;
      
    for (int i = 0; i < numParticles; i++)
    {
      h_x[i] = x[i];
      h_y[i] = y[i];
      h_z[i] = z[i];
      
      h_mass[i] = mass[i];
    }

    auto T_MovingData_6 = std::chrono::high_resolution_clock::now();
    cudaMemcpy(dev_x, h_x, numNodes * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y, h_y, numNodes * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_z, h_z, numNodes * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mass, h_mass, numNodes * sizeof(float), cudaMemcpyHostToDevice);
    auto end_MovingData_6 = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData_6 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData_6 - T_MovingData_6);
    cout << "T_MovingData_6 = " << elapsed_MovingData_6.count() * 1e-9 << endl;
    
    auto T_bounding_box = std::chrono::high_resolution_clock::now();
    compute_bounding_box_kernel<<< 1, blockSize_bh >>>(d_mutex, dev_x, dev_y, dev_z, d_left, d_right, d_bottom, d_top, d_front, d_back,
                                                                 numParticles);
    cudaDeviceSynchronize();
    auto end_bounding_box = std::chrono::high_resolution_clock::now();
    auto elapsed_bounding_box = std::chrono::duration_cast<std::chrono::nanoseconds>(end_bounding_box - T_bounding_box);
    cout << "T_bounding_box = " << elapsed_bounding_box.count() * 1e-9 << endl;

    T_build_tree_kernel = std::chrono::high_resolution_clock::now();
    build_tree_kernel<<< 1, 256 >>>(dev_x, dev_y, dev_z, dev_mass, d_count, d_start, d_child, d_index, d_left, d_right, d_bottom, d_top, d_front, d_back,
                                    numParticles, numNodes);
    cudaDeviceSynchronize();  
    end_build_tree_kernel = std::chrono::high_resolution_clock::now();
    elapsed_build_tree_kernel = std::chrono::duration_cast<std::chrono::nanoseconds>(end_build_tree_kernel - T_build_tree_kernel);
    cout << "Tree construction time = " << elapsed_build_tree_kernel.count() * 1e-9 << endl;
        
    centre_of_mass_kernel<<<1, blockSize_bh>>>(dev_x, dev_y, dev_z, dev_mass, d_index, numParticles);
    cudaDeviceSynchronize();  
    
    auto T_sort = std::chrono::high_resolution_clock::now();
    sort_kernel<<< 1, 256 >>>(d_count, d_start, d_sorted, d_child, d_index, numParticles);
    cudaDeviceSynchronize();
    auto end_sort = std::chrono::high_resolution_clock::now();
    auto elapsed_sort = std::chrono::duration_cast<std::chrono::nanoseconds>(end_sort - T_sort);
    cout << "T_sort = " << elapsed_sort.count() * 1e-9 << endl; 
    
    T_Force = std::chrono::high_resolution_clock::now();
    compute_forces_kernel<<< gridSize_bh, blockSize_bh >>>(dev_x, dev_y, dev_z, dev_ax, dev_ay, dev_az, dev_mass, d_eps, d_sorted, d_child,
                                                           d_left, d_right, d_bottom, d_top, d_front, d_back, numParticles);
    cudaDeviceSynchronize();
    end_Force = std::chrono::high_resolution_clock::now();
    elapsed_Force = std::chrono::duration_cast<std::chrono::nanoseconds>(end_Force - T_Force);
    cout << "T_Force = " << elapsed_Force.count() * 1e-9 << endl;
    
    
    auto T_MovingData_7 = std::chrono::high_resolution_clock::now();
    cudaMemcpy(h_ax, dev_ax, numNodes * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_ay, dev_ay, numNodes * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_az, dev_az, numNodes * sizeof(float), cudaMemcpyDeviceToHost);
    auto end_MovingData_7 = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData_7 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData_7 - T_MovingData_7);
    cout << "T_MovingData_7 = " << elapsed_MovingData_7.count() * 1e-9 << endl;

    for (int i = 0; i < numParticles; i++)
    {
      accx[i] = h_ax[i];
      accy[i] = h_ay[i];
      accz[i] = h_az[i];
    }

    auto T_MovingData_8 = std::chrono::high_resolution_clock::now();
    cudaMemcpy(d_accx, accx, numParticles * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_accy, accy, numParticles * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_accz, accz, numParticles * sizeof(float), cudaMemcpyHostToDevice);
    auto end_MovingData_8 = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData_8 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData_8 - T_MovingData_8);
    cout << "T_MovingData_8 = " << elapsed_MovingData_8.count() * 1e-9 << endl;

    delete[] h_mass;
    delete[] h_x;
    delete[] h_y;
    delete[] h_z;
    delete[] h_ax;
    delete[] h_ay;
    delete[] h_az;
    delete[] h_child;
    delete[] h_start;
    delete[] h_sorted;
    delete[] h_count;
    
    cudaFree(d_left);
    cudaFree(d_right);
    cudaFree(d_bottom);
    cudaFree(d_top);
    cudaFree(d_front);
    cudaFree(d_back);
    
    cudaFree(dev_mass);
    cudaFree(dev_x);
    cudaFree(dev_y);
    cudaFree(dev_z);
    
    cudaFree(dev_ax);
    cudaFree(dev_ay);
    cudaFree(dev_az);
    
    cudaFree(d_index);
    cudaFree(d_child);
    cudaFree(d_start);
    cudaFree(d_sorted);
    cudaFree(d_count);
    cudaFree(d_mutex);

    auto end_acc_g_Full_process = std::chrono::high_resolution_clock::now();
    auto elapsed_acc_g_Full_process = std::chrono::duration_cast<std::chrono::nanoseconds>(end_acc_g_Full_process - T_acc_g_Full_process);
    cout << "T_acc_g_Full_process = " << elapsed_acc_g_Full_process.count() * 1e-9 << endl;

    //****************** acc_tot **************************
    auto T_acc_tot = std::chrono::high_resolution_clock::now();
    acc_g_sph<<<gridSize, blockSize>>>(d_Typ, d_accx_tot, d_accy_tot, d_accz_tot,
                                       d_accx, d_accy, d_accz,
                                       d_accx_sph, d_accy_sph, d_accz_sph,
                                       N);
    cudaDeviceSynchronize();
    auto end_acc_tot = std::chrono::high_resolution_clock::now();
    auto elapsed_acc_tot = std::chrono::duration_cast<std::chrono::nanoseconds>(end_acc_tot - T_acc_tot);
    cout << "T_acc_tot = " << elapsed_acc_tot.count() * 1e-9 << endl;

    
    //******* Isothermal Gravity (Richings et al - 2018) ********
    galaxy_isothermal_potential<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_accx_tot,
                                                         d_accy_tot, d_accz_tot, sigma, G, N);
    cudaDeviceSynchronize();
    

    //****************** velocity evolution *******************
    v_evolve<<<gridSize, blockSize>>>(d_Typ, d_vx, d_vy, d_vz, d_accx_tot, d_accy_tot,
                                      d_accz_tot, dt, N);
    cudaDeviceSynchronize();

    //******************** get_dUX (du_dt) *********************
    auto T_dU = std::chrono::high_resolution_clock::now();
    get_dU_ngb<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                      d_divV, d_curlV, d_mass, d_P, d_dudt,
                                      d_ngb, MAX_ngb, visc_alpha, N);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) 
    {
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    cout << "get_dU ----> Successful!." << endl;
    auto end_dU = std::chrono::high_resolution_clock::now();
    auto elapsed_dU = std::chrono::duration_cast<std::chrono::nanoseconds>(end_dU - T_dU);
    cout << "T_dU = " << elapsed_dU.count() * 1e-9 << endl;

    //******************** u evolution *********************
    
    u_updater<<<gridSize, blockSize>>>(d_Typ, d_u, d_dudt, d_utprevious, dt, N);
    cudaDeviceSynchronize();
    
    
    //****************** Heating & Cooling - CLOUDY ********************
    auto T_cool = std::chrono::high_resolution_clock::now();
    doCloudyHCooling<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_rho, d_u, d_NH, d_nHGrid, d_TGrid, d_rGrid, d_NHGrid, d_uarr, d_muarr,
                                              d_Heat, d_Cool, N_nH, N_T, N_r, N_NH, gamma, kB, mH, nHLowBound, nHUpBound, TLowBound, TUpBound,
                                              rLowBound, rUpBound, NHLowBound, NHUpBound, UnitDensity_in_cgs, Unit_u_in_cgs, unitTime_in_s, dt, N);
    cudaDeviceSynchronize();
    auto end_cool = std::chrono::high_resolution_clock::now();
    auto elapsed_cool = std::chrono::duration_cast<std::chrono::nanoseconds>(end_cool - T_cool);
    cout << "T_cool = " << elapsed_cool.count() * 1e-9 << endl;
    
    
    //------ Cooling Test used for debugging -------------
    /*
    for (int i = 1110; i < 1115; i++)
    {
      cout << "u Before[" << i << "] = " << u[i] << endl;
    }
    
    cudaMemcpy(x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(y, d_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(z, d_z, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(rho, d_rho, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(NH, d_NH, N * sizeof(float), cudaMemcpyDeviceToHost);
    
    int ii = 1110;
    
    float uu = u[ii] * Unit_u_in_cgs;
    float rr = sqrt(x[ii]*x[ii] + y[ii]*y[ii] + z[ii]*z[ii]);
    float NHH = NH[ii];
    float nHH = rho[ii] * UnitDensity_in_cgs * XH / mH;
    
    printf("u, rr, NH, nH = %E, %f, %f, %f\n\n", uu, rr, NHH, nHH);
    
    cudaMemcpy(u, d_u, N * sizeof(float), cudaMemcpyDeviceToHost);
    
    for (int i = 1110; i < 1115; i++)
    {
      cout << "u After[" << i << "] = " << u[i] << endl;
    }
    
    cout << "Here we have an exit !!!" << endl;
    exit(0);
    */


    //-------------------------------------------------
    cudaMemcpy(rho, d_rho, N * sizeof(float), cudaMemcpyDeviceToHost);
    for (int i = 0; i < 5; i++)
    {
      cout << "AAA = " << rho[i] << endl;
    }

    auto T_SaveFile = std::chrono::high_resolution_clock::now();
    //------------ SAVING SNAP-SHOTS ------------
    cudaMemcpy(h, d_h, N * sizeof(float), cudaMemcpyDeviceToHost); // Moved outside so that it can be used by nSplit calculator in ach time step.
    if (!(counter % 200))
    //if (counter > -1)
    {
      cudaMemcpy(Typ, d_Typ, N * sizeof(int), cudaMemcpyDeviceToHost);

      cudaMemcpy(x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(y, d_y, N * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(z, d_z, N * sizeof(float), cudaMemcpyDeviceToHost);

      cudaMemcpy(vx, d_vx, N * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(vy, d_vy, N * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(vz, d_vz, N * sizeof(float), cudaMemcpyDeviceToHost);

      cudaMemcpy(rho, d_rho, N * sizeof(float), cudaMemcpyDeviceToHost);
      
      cudaMemcpy(NH, d_NH, N * sizeof(float), cudaMemcpyDeviceToHost);

      cudaMemcpy(u, d_u, N * sizeof(float), cudaMemcpyDeviceToHost);
      
      cudaMemcpy(mass, d_mass, N * sizeof(float), cudaMemcpyDeviceToHost);

      // Specify the output file name
      std::string filename = "./Outputs/G-" + to_string(t * 10) + ".bin";
      // Save the arrays to binary format
      saveArraysToBinary(filename, x, y, z, vx, vy, vz, rho, h, u, mass, NH, Typ, N);
    }
    auto end_SaveFile = std::chrono::high_resolution_clock::now();
    auto elapsed_SaveFile = std::chrono::duration_cast<std::chrono::nanoseconds>(end_SaveFile - T_SaveFile);
    cout << "T_SaveFile = " << elapsed_SaveFile.count() * 1e-9 << endl;


    //******************************************************
    //************* Updating Time-step dt ******************
    //******************************************************
    /*
    dt_array_indiv_dt<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z,
                                               d_vx, d_vy, d_vz,
                                               d_accx, d_accy, d_accz,
                                               d_accx_tot, d_accy_tot, d_accz_tot,
                                               d_h, d_csnd, d_dt_particles,
                                               d_abs_acc_g, d_abs_acc_tot,
                                               d_divV, d_dh_dt, C_CFL,
                                               visc_alpha, d_eps, N);
    cudaDeviceSynchronize();

    cudaMemcpy(dt_particles, d_dt_particles, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Typ, d_Typ, N * sizeof(int), cudaMemcpyDeviceToHost);
    */

    t += dt;

    // dt = min_finder(Typ, dt_particles, N);

    //***********************************************************
    //*************** Outflow particle injection ****************
    //***********************************************************
    // Generate a seed using the high resolution clock
    auto now = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    unsigned long long seed = counter; //static_cast<unsigned long long>(nanos); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //------------
    outflow_injector2<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z,
                                               d_vx, d_vy, d_vz,
                                               d_h, d_eps, d_mass,
                                               Nngb_f, d_Nngb_previous,
                                               d_u, M_dot_in, v_in,
                                               m_sph_high_res, u_for_10K_Temp,
                                               d_leftover_mass, dt, ndx_BH, N,
                                               seed);
    cudaDeviceSynchronize();

    if (!(counter % 1))
    {
      cout << "Adopted dt = " << dt << endl;
      cout << "current t = " << t << endl;
      cout << "*****************************" << endl;
      cout << endl;
    }

    delete[] ngb;
    cudaFree(d_ngb);
    delete[] groupedIndex;
    delete[] countx;
    delete[] offSet;
    
    cudaFree(d_groupedIndex);
    cudaFree(d_countx);
    cudaFree(d_offSet);

    counter++;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    cout << "Elapsed time = " << elapsed.count() * 1e-9 << endl;
    cout << endl;
    
  }

  delete[] Typ;
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] mass;
  delete[] h;
  delete[] rho;
  delete[] accx;
  delete[] accy;
  delete[] accz;
  delete[] eps;
  delete[] P;
  delete[] csnd;
  delete[] divV;
  delete[] curlV;
  delete[] accx_sph;
  delete[] accy_sph;
  delete[] accz_sph;
  delete[] accx_tot;
  delete[] accy_tot;
  delete[] accz_tot;
  delete[] abs_acc_g;
  delete[] abs_acc_tot;
  delete[] v_sig;
  delete[] dh_dt;
  delete[] u;
  delete[] dudt;
  delete[] utprevious;

  cudaFree(d_Typ);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_vx);
  cudaFree(d_vy);
  cudaFree(d_vz);
  cudaFree(d_mass);
  cudaFree(d_h);
  cudaFree(d_rho);
  cudaFree(d_accx);
  cudaFree(d_accy);
  cudaFree(d_accz);
  cudaFree(d_P);
  cudaFree(d_csnd);
  cudaFree(d_divV);
  cudaFree(d_curlV);
  cudaFree(d_accx_sph);
  cudaFree(d_accy_sph);
  cudaFree(d_accz_sph);
  cudaFree(d_accx_tot);
  cudaFree(d_accy_tot);
  cudaFree(d_accz_tot);
  cudaFree(d_abs_acc_g);
  cudaFree(d_abs_acc_tot);
  cudaFree(d_v_sig);
  cudaFree(d_dh_dt);
  cudaFree(d_u);
  cudaFree(d_dudt);
  cudaFree(d_utprevious);
}
