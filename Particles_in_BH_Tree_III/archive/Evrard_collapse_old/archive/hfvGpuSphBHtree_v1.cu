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
#include "hfvCLibs_v7.h"
#include "bh_tree_iteration_v2.h"
#include <cstdlib> // This is ONLY used for the "exit(0)" function !!

// Added the isothermal gravitational field acceleration. (24 May 2023).
// Added the reading of the params.txt file and updated the IC reading file section and function. (22 May 2023).

using namespace std;

int main()
{

  float dt = 1e-4; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This is only the first time step !!

  const float Nngb_f = 64.0f; // used in smoothing func.
  const int Nngb = 64;
  const int Ndown = Nngb - 5;
  const int Nup = Nngb + 5;
  const float coeff = 0.005f; // used for smoothing length.
  
  const float kpc_in_cm = 3.086e21;
  
  //*******************************************************************
  //******************* Reading Cooling File **************************
  //*******************************************************************
  ifstream infile("coolHeatGridZ.bin", ios::binary);

  if (!infile) {
    cerr << "Failed to open coolHeatGrid.bin file." << endl;
    return 1;
  }

  // Read the sizes
  int N_kpc, N_T, N_nH, N_Z, N_M, N_Time;
  infile.read(reinterpret_cast<char*>(&N_kpc), sizeof(int));
  infile.read(reinterpret_cast<char*>(&N_nH), sizeof(int));
  infile.read(reinterpret_cast<char*>(&N_Z), sizeof(int));
  infile.read(reinterpret_cast<char*>(&N_T), sizeof(int));
  infile.read(reinterpret_cast<char*>(&N_M), sizeof(int));
  infile.read(reinterpret_cast<char*>(&N_Time), sizeof(int));

  // Allocate and read the densities, temperatures, metallicities, and timeArr arrays
  vector<float> kpcArr(N_kpc);     // float
  vector<float> densities(N_nH);     // float
  vector<float> metallicities(N_Z);  // float
  vector<float> temperatures(N_T);   // float
  vector<float> timeArr(N_Time);     // float

  infile.read(reinterpret_cast<char*>(kpcArr.data()), N_kpc * sizeof(float));
  infile.read(reinterpret_cast<char*>(densities.data()), N_nH * sizeof(float));
  infile.read(reinterpret_cast<char*>(metallicities.data()), N_Z * sizeof(float));
  infile.read(reinterpret_cast<char*>(temperatures.data()), N_T * sizeof(float));
  infile.read(reinterpret_cast<char*>(timeArr.data()), N_Time * sizeof(float));

  // Allocate and read the flattened res and muArr array
  int N_HCool = N_kpc * N_T * N_nH * N_Z * N_Time;
  vector<float> res_flattened(N_HCool);  // float
  vector<float> muArr(N_HCool);  // float
  
  int N_metalz = N_kpc * N_T * N_nH * N_Z * N_M * N_Time;
  vector<float> metalzArr(N_metalz);  // float

  infile.read(reinterpret_cast<char*>(res_flattened.data()), N_HCool * sizeof(float));
  infile.read(reinterpret_cast<char*>(muArr.data()), N_HCool * sizeof(float)); // Note that muA and res_flattedned have the same structure!!
  
  infile.read(reinterpret_cast<char*>(metalzArr.data()), N_metalz * sizeof(float));

  infile.close();
  
  
//------------- Just for testing ---------
  int jjj = 2; // kpc
  int i = 20;   // T
  int j = 70;  // nH
  int k = 2;   // Z
  int l = 10;   // time
  
  int indx = jjj * (N_T * N_nH * N_Z * N_Time) + i * (N_nH * N_Z * N_Time) + j * (N_Z * N_Time) + k * N_Time + l;
  
  int ii_HI  = 0;
  int ii_HII = 1;
  int indx_HI  = jjj * (N_T * N_nH * N_Z * N_M * N_Time) + i * (N_nH * N_Z * N_M * N_Time) + j * (N_Z * N_M * N_Time) + k * (N_M * N_Time) + ii_HI * (N_Time) + l;
  int indx_HII = jjj * (N_T * N_nH * N_Z * N_M * N_Time) + i * (N_nH * N_Z * N_M * N_Time) + j * (N_Z * N_M * N_Time) + k * (N_M * N_Time) + ii_HII * (N_Time) + l;
  
  cout << "u = " << res_flattened[indx] << endl;
  cout << "mu = " << muArr[indx] << endl;
  cout << "HI fraction = " << metalzArr[indx_HI] << endl;
  cout << "HII fraction = " << metalzArr[indx_HII] << endl;
//--------------------------------------------
  
  float *Temp, *d_Temp, *nH, *d_nH, *Z, *d_Z, *Time, *d_Time, *HCool, *d_HCool, *muA, *d_muA, *kpc, *d_kpc, *metalz, *d_metalz;
  float *U, *d_U; // They will be used only inside the hcooling function as I could not define them inside the hcooling function (GPU climitations!!)
  
  kpc = new float[N_kpc];
  Temp = new float[N_T];
  nH = new float[N_nH];
  Z = new float[N_Z];
  Time = new float[N_Time];
  HCool = new float[N_HCool];
  muA = new float[N_HCool];
  U = new float[N_T];
  metalz = new float[N_metalz];
  
  cudaMalloc(&d_kpc, N_kpc * sizeof(float));
  cudaMalloc(&d_Temp, N_T * sizeof(float));
  cudaMalloc(&d_nH, N_nH * sizeof(float));
  cudaMalloc(&d_Z, N_Z * sizeof(float));
  cudaMalloc(&d_Time, N_Time * sizeof(float));
  cudaMalloc(&d_HCool, N_HCool * sizeof(float));
  cudaMalloc(&d_muA, N_HCool * sizeof(float));
  cudaMalloc(&d_U, N_T * sizeof(float));
  cudaMalloc(&d_metalz, N_metalz * sizeof(float));
  
  for (int i = 0; i < N_kpc; i++)
  {
    kpc[i] = kpcArr[i];
  }
  
  for (int i = 0; i < N_T; i++)
  {
    Temp[i] = temperatures[i];
    U[i] = 0.0f;
  }
  
  for (int i = 0; i < N_nH; i++)
  {
    nH[i] = densities[i];
  }
  
  for (int i = 0; i < N_Z; i++)
  {
    Z[i] = metallicities[i];
  }
  
  for (int i = 0; i < N_Time; i++)
  {
    Time[i] = timeArr[i];
  }
  
  for (int i = 0; i < N_HCool; i++)
  {
    HCool[i] = res_flattened[i];
  }
  
  for (int i = 0; i < N_HCool; i++)
  {
    muA[i] = muArr[i];
  }

  for (int i = 0; i < N_metalz; i++)
  {
    metalz[i] = metalzArr[i];
  }

  // Copy from Host to Device
  cudaMemcpy(d_kpc, kpc, N_kpc * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Temp, Temp, N_T * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nH, nH, N_nH * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Z, Z, N_Z * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Time, Time, N_Time * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_HCool, HCool, N_HCool * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_muA, muA, N_HCool * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_U, U, N_T * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_metalz, metalz, N_metalz * sizeof(float), cudaMemcpyHostToDevice);

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

  float gammah = 5.0f / 3.0f;
  float GAMMA_MINUS1 = gammah - 1.0f;
  
  int N_ionFrac = N * N_M; // We have N_M species for each particle (N = total number of particles)
  
  float *ionFrac, *d_ionFrac;

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
  
  ionFrac = new float[N_ionFrac];

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
  
  cudaMalloc(&d_ionFrac, N_ionFrac * sizeof(float));

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
  
  for (int i = 0; i < N_ionFrac; i++)
  {
    ionFrac[i] = 0.0;
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
  
  cudaMemcpy(d_ionFrac, ionFrac, N_ionFrac * sizeof(float), cudaMemcpyHostToDevice);
  
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

  const float visc_alpha = 1.0f;

  float t;

  t = 0.0f;

  float tEnd = 1.0f;
  float Nt = ceil(tEnd / dt) + 1;

  float Zmetal = 0.1; // ==> [Z/H] = -1.

  //-----------------------------------------------
  //-------------- Smoothing Length ---------------
  //-----------------------------------------------
  smoothing_h<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_h,
                                       N, Ndown, Nup, coeff,
                                       Nngb_f, d_Nngb_previous, d_divV, dt);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //----------------- getDensity ------------------
  //-----------------------------------------------
  getDensity<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_mass,
                                      d_rho, d_h, N);
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
  div_curlVel<<<gridSize, blockSize>>>(d_Typ, d_divV, d_curlV, d_x, d_y, d_z, d_vx, d_vy, d_vz,
                                       d_rho, d_mass, d_h, N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //------------------ acc_sph --------------------
  //-----------------------------------------------
  acc_sph<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                   d_divV, d_curlV, d_mass, d_P, d_accx_sph, d_accy_sph,
                                   d_accz_sph, visc_alpha, N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //------------------ getAcc_g -------------------
  //-----------------------------------------------
  /*
  acc_g<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_eps, d_accx, d_accy, d_accz,
                                 d_mass, G, N);
  cudaDeviceSynchronize();
  */
  
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
  numNodes = 8 * numParticles + 10000000;

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
  
  //exit(0);

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

  get_dU<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                  d_divV, d_curlV, d_mass, d_P, d_dudt,
                                  visc_alpha, N);
  cudaDeviceSynchronize();

  //-----------------------------------------------
  //---------------- u evolution ------------------
  //-----------------------------------------------

  u_updater<<<gridSize, blockSize>>>(d_Typ, d_u, d_dudt, d_utprevious, dt, N);
  cudaDeviceSynchronize();

  //const float C_CFL = 0.25;

  //float h_min, h_max, h_mean;

  float leftover_mass = 0.0f;
  float *d_leftover_mass;
  cudaMalloc((void **)&d_leftover_mass, sizeof(float));
  cudaMemcpy(d_leftover_mass, &leftover_mass, sizeof(float), cudaMemcpyHostToDevice);

  // **************************************************************
  // *********************** MAIN LOOP ****************************
  // **************************************************************

  int counter = 0; // This is used to save fewer output files, e.g. 1 snap-shot per 2 time-step!

  while (t < tEnd)
  {

    auto begin = std::chrono::high_resolution_clock::now();

    //****************** velocity evolution *******************
    v_evolve<<<gridSize, blockSize>>>(d_Typ, d_vx, d_vy, d_vz, d_accx_tot, d_accy_tot,
                                      d_accz_tot, dt, N);
    cudaDeviceSynchronize();

    //****************** position evolution (BH fixed at [0, 0, 0]) *******************

    r_evolveT<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, dt, N);
    cudaDeviceSynchronize();

    //****************** Smoothing Length *********************
    auto T_hh = std::chrono::high_resolution_clock::now();
    smoothing_h<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_h,
                                         N, Ndown, Nup, coeff,
                                         Nngb_f, d_Nngb_previous, d_divV, dt);
    cudaDeviceSynchronize();
    auto end_hh = std::chrono::high_resolution_clock::now();
    auto elapsed_hh = std::chrono::duration_cast<std::chrono::nanoseconds>(end_hh - T_hh);
    cout << "T_h = " << elapsed_hh.count() * 1e-9 << endl;

    //****************** Set eps of Gas equal to h ******************
    set_eps_of_gas_to_h<<<gridSize, blockSize>>>(d_Typ, d_eps, d_h, N);
    cudaDeviceSynchronize();

    //****************** getDensity ***********************
    auto T_density = std::chrono::high_resolution_clock::now();
    getDensity<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_mass,
                                        d_rho, d_h, N);
    cudaDeviceSynchronize();
    auto end_density = std::chrono::high_resolution_clock::now();
    auto elapsed_density = std::chrono::duration_cast<std::chrono::nanoseconds>(end_density - T_density);
    cout << "T_density = " << elapsed_density.count() * 1e-9 << endl;
    
    
    //****************** getPressure **********************
    getPressure_Adiabatic<<<gridSize, blockSize>>>(d_Typ, d_P, d_rho, d_u, gammah, N);
    cudaDeviceSynchronize();

    //****************** getCsound ************************
    getCsound_Adiabatic<<<gridSize, blockSize>>>(d_Typ, d_csnd, d_u, gammah, N);
    cudaDeviceSynchronize();

    //****************** div_curlVX ************************
    auto T_divCurl = std::chrono::high_resolution_clock::now();
    div_curlVel<<<gridSize, blockSize>>>(d_Typ, d_divV, d_curlV, d_x, d_y, d_z, d_vx, d_vy, d_vz,
                                         d_rho, d_mass, d_h, N);
    cudaDeviceSynchronize();
    auto end_divCurl = std::chrono::high_resolution_clock::now();
    auto elapsed_divCurl = std::chrono::duration_cast<std::chrono::nanoseconds>(end_divCurl - T_divCurl);
    cout << "T_divCurl = " << elapsed_divCurl.count() * 1e-9 << endl;

    //****************** acc_sphX **************************
    auto T_acc_sph = std::chrono::high_resolution_clock::now();
    acc_sph<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                     d_divV, d_curlV, d_mass, d_P, d_accx_sph, d_accy_sph,
                                     d_accz_sph, visc_alpha, N);
    cudaDeviceSynchronize();
    auto end_acc_sph = std::chrono::high_resolution_clock::now();
    auto elapsed_acc_sph = std::chrono::duration_cast<std::chrono::nanoseconds>(end_acc_sph - T_acc_sph);
    cout << "T_acc_sph = " << elapsed_acc_sph.count() * 1e-9 << endl;
    
    
    
    //****************** getAcc_gX *************************
    /*
    auto T_acc_g = std::chrono::high_resolution_clock::now();
    acc_g<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_eps, d_accx, d_accy, d_accz,
                                   d_mass, G, N);
    cudaDeviceSynchronize();
    auto end_acc_g = std::chrono::high_resolution_clock::now();
    auto elapsed_acc_g = std::chrono::duration_cast<std::chrono::nanoseconds>(end_acc_g - T_acc_g);
    cout << "T_acc_g = " << elapsed_acc_g.count() * 1e-9 << endl;
    */

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
    numNodes = 8 * numParticles + 10000000;

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
   
    auto T_reset_kernel = std::chrono::high_resolution_clock::now();
    reset_arrays_kernel<<< 1, blockSize_bh >>>(d_mutex, dev_x, dev_y, dev_z, dev_mass, d_count, d_start, d_sorted, d_child, d_index,
                                                         d_left, d_right, d_bottom, d_top, d_front, d_back, numParticles, numNodes);
    cudaDeviceSynchronize();
    auto end_reset_kernel = std::chrono::high_resolution_clock::now();
    auto elapsed_reset_kernel = std::chrono::duration_cast<std::chrono::nanoseconds>(end_reset_kernel - T_reset_kernel);
    cout << "T_reset_kernel = " << elapsed_reset_kernel.count() * 1e-9 << endl;
    
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
    
    auto T_bounding_box = std::chrono::high_resolution_clock::now();
    compute_bounding_box_kernel<<< 1, blockSize_bh >>>(d_mutex, dev_x, dev_y, dev_z, d_left, d_right, d_bottom, d_top, d_front, d_back,
                                                                 numParticles);
    cudaDeviceSynchronize();
    auto end_bounding_box = std::chrono::high_resolution_clock::now();
    auto elapsed_bounding_box = std::chrono::duration_cast<std::chrono::nanoseconds>(end_bounding_box - T_bounding_box);
    cout << "T_bounding_box = " << elapsed_bounding_box.count() * 1e-9 << endl;
    
    
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
    
    
    auto T_MovingData = std::chrono::high_resolution_clock::now();
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
    auto end_MovingData = std::chrono::high_resolution_clock::now();
    auto elapsed_MovingData = std::chrono::duration_cast<std::chrono::nanoseconds>(end_MovingData - T_MovingData);
    cout << "T_MovingData = " << elapsed_MovingData.count() * 1e-9 << endl;


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

    /*
    //******* Isothermal Gravity (Richings et al - 2018) ********
    galaxy_isothermal_potential<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_accx_tot,
                                                         d_accy_tot, d_accz_tot, sigma, G, N);
    cudaDeviceSynchronize();
    */

    //****************** velocity evolution *******************
    v_evolve<<<gridSize, blockSize>>>(d_Typ, d_vx, d_vy, d_vz, d_accx_tot, d_accy_tot,
                                      d_accz_tot, dt, N);
    cudaDeviceSynchronize();

    //******************** get_dUX (du_dt) *********************
    auto T_dU = std::chrono::high_resolution_clock::now();
    get_dU<<<gridSize, blockSize>>>(d_Typ, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_h, d_csnd, d_rho,
                                    d_divV, d_curlV, d_mass, d_P, d_dudt,
                                    visc_alpha, N);
    cudaDeviceSynchronize();
    auto end_dU = std::chrono::high_resolution_clock::now();
    auto elapsed_dU = std::chrono::duration_cast<std::chrono::nanoseconds>(end_dU - T_dU);
    cout << "T_dU = " << elapsed_dU.count() * 1e-9 << endl;

    //******************** u evolution *********************
    
    u_updater<<<gridSize, blockSize>>>(d_Typ, d_u, d_dudt, d_utprevious, dt, N);
    cudaDeviceSynchronize();
    
    
    /*
    //****************** Heating & Cooling ********************
    
    auto T_cool = std::chrono::high_resolution_clock::now();
    hcoolingx<<<gridSize, blockSize>>>(d_Typ, d_u, d_U, d_rho, d_metalz, Zmetal, dt, // Zmetal is the gass metallicity assumed.
                                      d_nH, d_Z, d_HCool, d_ionFrac, d_Time, d_x, d_y, d_z,
                                      d_muA, d_Temp, d_kpc, UnitDensity_in_cgs, unitTime_in_s, 
                                      Unit_u_in_cgs, unitLength_in_cm, kpc_in_cm, GAMMA_MINUS1,
                                      N_kpc, N_nH, N_Z, N_T, N_M, N_Time, N);
    cudaDeviceSynchronize();
    auto end_cool = std::chrono::high_resolution_clock::now();
    auto elapsed_cool = std::chrono::duration_cast<std::chrono::nanoseconds>(end_cool - T_cool);
    cout << "T_cool = " << elapsed_cool.count() * 1e-9 << endl;
    */
    
    //-------------------------------------------------

    cudaMemcpy(rho, d_rho, N * sizeof(float), cudaMemcpyDeviceToHost);
    for (int i = 0; i < 5; i++)
    {
      cout << "AAA = " << rho[i] << endl;
    }

    //------------ SAVING SNAP-SHOTS ------------
    if (!(counter % 100))
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
      cudaMemcpy(h, d_h, N * sizeof(float), cudaMemcpyDeviceToHost);

      cudaMemcpy(u, d_u, N * sizeof(float), cudaMemcpyDeviceToHost);
      
      cudaMemcpy(mass, d_mass, N * sizeof(float), cudaMemcpyDeviceToHost);

      cudaMemcpy(ionFrac, d_ionFrac, N_ionFrac * sizeof(float), cudaMemcpyDeviceToHost);

      // Specify the output file name
      std::string filename = "./Outputs/G-" + to_string(t * 1) + ".bin";
      // Save the arrays to binary format
      saveArraysToBinary(filename, x, y, z, vx, vy, vz, rho, h, u, mass, ionFrac, Typ, N, N_ionFrac);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    cout << "Elapsed time = " << elapsed.count() * 1e-9 << endl;
    cout << endl;

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

    /*
    // Generate a seed using the high resolution clock
    auto now = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    unsigned long long seed = static_cast<unsigned long long>(nanos);
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
    */

    if (!(counter % 1))
    {
      cout << "Adopted dt = " << dt << endl;
      cout << "current t = " << t << endl;
      cout << "*****************************" << endl;
      cout << endl;
    }

    counter++;
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
