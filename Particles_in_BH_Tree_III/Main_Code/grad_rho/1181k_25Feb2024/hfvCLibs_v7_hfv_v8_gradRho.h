#ifndef HFVCPPLIBS_H
#define HFVCPPLIBS_H

#include <curand.h>
#include <curand_kernel.h>

// ngbDebug is removed as it was useless after I defined Lim in ngbFinder (5 Jan 2024).
// I set a lower limit for u in the hcooling function corresponding to ~ 100 K (4 Jan 2024).
// I added a condition to prevent u from becoming negative after cooling/adiabatic expansion (3 Jan 2024).
// I debugged the fast cooling problem and implemented the new hcool function (27 Dec 2023).
// I used my method of grid construction (ngbDB_v2) to speed up the code (1 Nov 2023).
// smoothing_h function updated (16 Oct 2023).
// Updating the "h" assignment for outflow particles to improve the smoothin_h calculation speed (13 Oct 2023).
// We Re-set the BH position back to (0, 0, 0). This was forgotten to be implemented !! (13 Oct 2023).
// phi updated in outflow injection function (4 Oct 2023).
// Heating, Cooling using CHIMES added. Also UnitDensity_in_cgs, Unit_u_in_cgs, unitTime_in_s added to params.txt (11 Aug 2023).
// Added the isothermal gravitational field acceleration. (24 May 2023).
// Added the reading of the params.txt file and updated the IC reading file section and function. (22 May 2023).

//******************************************
//********** Reading params.txt ************
//******************************************

void readParams(std::string &filename, int &N_tot, int &ndx_BH, float &G, float &L_AGN_code_unit,
                float &M_dot_in_code_unit, float &vin_in_code_unit,
                float &u_for_10K_Temp, float &m_sph_high_res, float &sigma,
                float &UnitDensity_in_cgs, float &Unit_u_in_cgs, float &unitTime_in_s,
                float &unitLength_in_cm)
{
  std::ifstream file("params.txt");
  if (file.is_open())
  {
    std::getline(file, filename); // Read filename string
    file >> N_tot;                // Read N_tot
    file >> ndx_BH;               // Read ndx_BH
    file >> G;                    // Read G
    file >> L_AGN_code_unit;      // Read L_AGN_code_unit
    file >> M_dot_in_code_unit;   // Read M_dot_in_code_unit
    file >> vin_in_code_unit;     // Read vin_in_code_unit
    file >> u_for_10K_Temp;       // Read u_for_10K_Temp
    file >> m_sph_high_res;       // Read m_sph_high_res
    file >> sigma;                // Read sigma
    file >> UnitDensity_in_cgs;                // Read UnitDensity_in_cgs
    file >> Unit_u_in_cgs;                // Read Unit_u_in_cgs
    file >> unitTime_in_s;                // Read unitTime_in_s
    file >> unitLength_in_cm;                // Read unitLength_in_cm
  }
  else
  {
    std::cout << "Unable to open params.txt file";
  }
  file.close();
}

//***************************************
//********** Reading IC file ************
//***************************************
std::tuple<std::vector<int>, 
           std::vector<float>, std::vector<float>, std::vector<float>, 
           std::vector<float>, std::vector<float>, std::vector<float>, 
           std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>>
readVectorsFromFile(const std::string &filename) 
{
    int N;
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Could not open the IC file for reading." << std::endl;
        exit(1); // or other error handling
    }

    // Read N_tot from the start of the file
    file.read(reinterpret_cast<char*>(&N), sizeof(int));

    // Create and resize vectors based on N
    std::vector<int> Typvec(N);
    std::vector<float> xvec(N);
    std::vector<float> yvec(N);
    std::vector<float> zvec(N);
    std::vector<float> vxvec(N);
    std::vector<float> vyvec(N);
    std::vector<float> vzvec(N);
    std::vector<float> uvec(N);
    std::vector<float> hvec(N);
    std::vector<float> epsvec(N);
    std::vector<float> massvec(N);

    file.read(reinterpret_cast<char*>(Typvec.data()), sizeof(int) * N);
    file.read(reinterpret_cast<char*>(xvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(yvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(zvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(vxvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(vyvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(vzvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(uvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(hvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(epsvec.data()), sizeof(float) * N);
    file.read(reinterpret_cast<char*>(massvec.data()), sizeof(float) * N);

    file.close();

    return std::make_tuple(Typvec, xvec, yvec, zvec, vxvec, vyvec, vzvec, uvec, hvec, epsvec, massvec);
}



//*************************************************************************
//*************** Function to save the OUTPUT Snap-Shots!! ****************
//*************************************************************************
void saveArraysToBinary(const std::string &filename, float *x, float *y, float *z, float *vx, float *vy, float *vz,
                        float *rho, float *h, float *u, float *mass, float *gradRhoNorm, int *Typ, int N)
{
  // Open the file in binary mode
  std::ofstream file(filename, std::ios::binary);

  // Check if the file was opened successfully
  if (!file)
  {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  // Write N and NG to the file
  file.write(reinterpret_cast<const char *>(&N), sizeof(int));

  // Write the arrays to the file
  file.write(reinterpret_cast<const char *>(Typ), N * sizeof(int));
  file.write(reinterpret_cast<const char *>(x), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(y), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(z), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(vx), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(vy), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(vz), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(rho), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(h), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(u), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(mass), N * sizeof(float));
  file.write(reinterpret_cast<const char *>(gradRhoNorm), N * sizeof(float));

  // Close the file
  file.close();
}



//*******************************
//********* max_finder **********
//*******************************
float max_finder(int *Typ, float *arr, int N)
{

  float max_val = 0.0;
  for (int i = 0; i < N; i++)
  {
    if (Typ[i] == 0)
    {
      if (arr[i] >= max_val)
      {
        max_val = arr[i];
      }
    }
  }
  return max_val;
}

//*******************************
//********* min_finder **********
//*******************************
float min_finder(int *Typ, float *arr, int N)
{

  float min_val = 1e22; // used large value for the start! I did not want to use arr[0] as it could be 0.0 itself!
  for (int i = 0; i < N; i++)
  {
    if (Typ[i] == 0)
    {
      if (arr[i] <= min_val)
      {
        min_val = arr[i];
      }
    }
  }
  return min_val;
}


//========================================
//================ ngbDB =================
//========================================
__global__ void ngbDB(int *Typ, float *x, float *y, float *z, float *h, int *ngb, int MAX_ngb, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float h_new = 2.10f * h[i];

    int k = 0;
    int j = 0;

    float dx, dy, dz;
    float rr;

    while (j < N && k < MAX_ngb)
    {
      if (Typ[j] == 0)
      {
        dx = x[j] - x[i];
        dy = y[j] - y[i];
        dz = z[j] - z[i];
        rr = sqrt(dx * dx + dy * dy + dz * dz);

        if (rr <= h_new)
        {
          ngb[i * MAX_ngb + k] = j;
          k++;
        }
      }
      j++;
    }
  }
}




//========================================
//================ ngbDB_v2 =================
//========================================
__global__ void ngbDB_v2(int *Typ, float *x, float *y, float *z, float *h, int *ngb, int MAX_ngb, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float coeff = 0.05f;
    float counter = 0.0f; // is used to increase the search radius (i.e. h_new) in case the number of neighbrs is less than 65.
    
    float h_new = 0.0f;

    int k = 0;
    int j = 0;

    float dx, dy, dz;
    float rr;

    while (k < 65)
    {
      k = 0;
      
      h_new = (2.0f + (counter + 1.0f) * coeff) * h[i];
      
      j = 0;
      
      while (j < N && k < MAX_ngb)
      {
        if (Typ[j] == 0)
        {
          dx = x[j] - x[i];
          dy = y[j] - y[i];
          dz = z[j] - z[i];
          rr = sqrt(dx * dx + dy * dy + dz * dz);

          if (rr <= h_new)
          {
            ngb[i * MAX_ngb + k] = j;
            k++;
          }
        }
        j++;
      }
      counter++;
    }
  }
}




//========================================
//========== Smoothing Length ============ Updated: 28 Jan 2023. h_new adopted from eq.31 in Gadget2 Paper
//========================================
__global__ void smoothing_h(int *Typ, float *x, float *y, float *z, float *h,
                            int N, int Ndown, int Nup, float coeff,
                            float Nngb_f, float *Nngb_previous, float *divV, float dt)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float h_new = 2.0f * (0.5f * h[i] * (1.0f + pow((Nngb_f / Nngb_previous[i]), 1.0f / 3.0f)) +
                          1.0f / 3.0f * h[i] * divV[i] * dt);
    
    int N_iter = 0;
    int k = 0;
    
    int k_pre = Nup; // Just a choice so that it does not activate that if condition in the first run!

    float dx, dy, dz;
    while ((k < Ndown) || (k > Nup))
    {

      k = 0;

      for (int j = 0; j < N; j++)
      {
        if (Typ[j] == 0)
        {
          dx = x[j] - x[i];
          dy = y[j] - y[i];
          dz = z[j] - z[i];
          float rr = sqrt(dx * dx + dy * dy + dz * dz);

          if (rr <= h_new)
          {
            k++;
          }
        }
        
        //-----------
        if (k > Nup) // To stop unnecessary search after it reaches k > Nup!
        {
          break;
        }
        //-----------
      }

      //-----------
      if (((k < Ndown) && (k_pre > Nup)) || ((k > Nup) && (k_pre < Ndown))) // To prevent oscillation outside Nup and Ndown values!!
      {
        coeff = coeff / 2.0f;
      }
      //-----------

      if (k < Ndown)
      {
        h_new = h_new + coeff * 2.0f * h[i];
      }
      
      if (k > Nup)
      {
        h_new = h_new - coeff * 2.0f * h[i];
      }
      
      k_pre = k;

      N_iter++;
      if (N_iter > 500)
      {
        break;
      }
    }
    Nngb_previous[i] = k;
    h[i] = 0.5 * h_new;
  }
}


//========================================
//========== Smoothing Length_ngb ============ Updated: 28 Jan 2023. h_new adopted from eq.31 in Gadget2 Paper
//========================================
__global__ void smoothing_h_ngb(int *Typ, float *x, float *y, float *z, float *h,
                                int N, int Ndown, int Nup, float coeff,
                                float Nngb_f, float *Nngb_previous, float *divV, int *ngb, int MAX_ngb, float dt)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float h_new = 2.0f * (0.5f * h[i] * (1.0f + pow((Nngb_f / Nngb_previous[i]), 1.0f / 3.0f)) +
                          1.0f / 3.0f * h[i] * divV[i] * dt);
    int N_iter = 0;
    int k = 0;
    int j = 0;
    
    int k_pre = Nup; // Just a choice so that it does not activate the 'if' condition in the first run!

    float dx, dy, dz;
    while ((k < Ndown) || (k > Nup))
    {

      k = 0;

      for (int jj = 0; jj < MAX_ngb; jj++)
      {

        j = ngb[i * MAX_ngb + jj];
        
        if (j != -1) // Note that the Typ[j] == 0 condition was already checked during ngb array construction!
        {
          dx = x[j] - x[i];
          dy = y[j] - y[i];
          dz = z[j] - z[i];
          float rr = sqrt(dx * dx + dy * dy + dz * dz);
          
          if (rr <= h_new)
          {
            k++;
          }
        }

        //-----------
        if (k > Nup) // To stop unnecessary search after it reaches k > Nup!
        {
          break;
        }
        //-----------
      }

      //-----------
      if (((k < Ndown) && (k_pre > Nup)) || ((k > Nup) && (k_pre < Ndown))) // To prevent oscillation outside Nup and Ndown values!!
      {
        coeff = coeff / 2.0f;
      }
      //-----------

      if (k < Ndown)
      {
        h_new = h_new + coeff * 2.0f * h[i];
      }

      if (k > Nup)
      {
        h_new = h_new - coeff * 2.0f * h[i];
      }
      
      k_pre = k;

      N_iter++;
      if (N_iter > 200)
      {
        printf("Max N_iter REACHED in smoothing_h_ngb for i, k, h_new, h[i], x, y, z =  !!!!!!!!!!!!!!!!!!!!!!!!!!! %d, %d, %f, %f, %f, %f, %f\n",
                                                          i, k, h_new, h[i], x[i], y[i], z[i]);
        h_new = 2.0 * h[i]; // When max N_iter reached we do not change the smoothing length for safety!!
        break;
      }
    }
    
    Nngb_previous[i] = k;
    h[i] = 0.5 * h_new;
  }
}

//==========================================
//========== Set eps of Gas to h ===========
//==========================================
__global__ void set_eps_of_gas_to_h(int *Typ, float *eps, float *h, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    eps[i] = h[i];
  }
}

//==========================================
//============ getDensity_ngb ==============
//==========================================
__global__ void getDensity_ngb(int *Typ, float *x, float *y, float *z, float *mass,
                           float *rho, float *h, int *ngb, int MAX_ngb, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float dx, dy, dz, rr, hij, sig, q, hij3;
    float WIij;
    float ss = 0.0f;
    
    int j = 0;

    for (int jj = 0; jj < MAX_ngb; jj++)
    {
    
      j = ngb[i * MAX_ngb + jj];
    
      if (j != -1) // Note that the Typ[j] == 0 condition was already checked during ngb array construction!
      {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        dz = z[i] - z[j];

        rr = sqrt(dx * dx + dy * dy + dz * dz);
        hij = 0.5f * (h[i] + h[j]);

        if (rr <= 2.0f * hij)
        {

          sig = 1.0 / M_PI;
          q = rr / hij;
          hij3 = hij * hij * hij;
          WIij = 0.0f;

          if (q <= 1.0)
          {
            WIij = sig / hij3 * (1.0f - (3.0f / 2.0f) * q * q + (3.0f / 4.0f) * q * q * q);
          }

          if ((q > 1.0f) && (q <= 2.0))
          {
            WIij = sig / hij3 * (1.0f / 4.0f) * (2.0f - q) * (2.0f - q) * (2.0f - q);
          }

          ss += mass[j] * WIij;
        }
      }
    }
    rho[i] = ss;
  }
}

//==============================================
//================= acc_g ======================
//==============================================
__global__ void acc_g(int *Typ, float *x, float *y, float *z, float *eps, float *accx,
                      float *accy, float *accz, float *mass, float G, int N)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && ((Typ[i] == 0) || (Typ[i] == 1)))
  {

    float dx, dy, dz, rr, inv_r3, epsij, q, q2, q3, q4, q5, q6, fk;
    float accxt = 0.0f, accyt = 0.0f, acczt = 0.0f;
    for (int j = 0; j < N; j++)
    {
      if ((Typ[j] == 0) || (Typ[j] == 1))
      {
        dx = x[j] - x[i];
        dy = y[j] - y[i];
        dz = z[j] - z[i];

        rr = sqrt(dx * dx + dy * dy + dz * dz);
        inv_r3 = 1.0f / (rr * rr * rr + 1e-5);
        epsij = 0.5f * (eps[i] + eps[j]);
        q = rr / epsij;
        q2 = q * q;
        q3 = q2 * q;
        q4 = q3 * q;
        q5 = q4 * q;
        q6 = q5 * q;

        if (q <= 1.0f)
        {
          fk = (1.0f / (epsij * epsij * epsij)) * ((4.0f / 3.0f) - (6.0f / 5.0f) * q2 + (1.0f / 2.0f) * q3);
        }

        if ((q > 1.0f) && (q <= 2.0f))
        {
          fk = inv_r3 * ((-1.0f / 15.0f) + (8.0f / 3.0f) * q3 - 3.0f * q4 + (6.0f / 5.0f) * q5 - (1.0f / 6.0f) * q6);
        }

        if (q > 2.0f)
        {
          fk = inv_r3;
        }

        accxt += G * fk * dx * mass[j];
        accyt += G * fk * dy * mass[j];
        acczt += G * fk * dz * mass[j];
      }
    }
    accx[i] = accxt;
    accy[i] = accyt;
    accz[i] = acczt;
  }
}

//===================================================
//============= getPressure (Adiabatic) =============
//===================================================
__global__ void getPressure_Adiabatic(int *Typ, float *P, float *rho, float *u, float gammah, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    P[i] = (gammah - 1.0f) * rho[i] * u[i];
  }
}

//===============================================
//============= getCsound (Adiabatic) ===========
//===============================================
__global__ void getCsound_Adiabatic(int *Typ, float *csnd, float *u, float gammah, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    csnd[i] = sqrt(gammah * (gammah - 1.0f) * u[i]);
  }
}



//=========================================================
//================== grad_rho_norm_ngb ====================
//=========================================================
__global__ void grad_rho_norm_ngb(int *Typ, float *x, float *y, float *z, float *rho, float *mass,
                                  int *ngb, int MAX_ngb, float *h, float *gradRhoNorm, int N)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float dx, dy, dz, rr, hij, q, hij5, sig;
    float nW = 0.0f;
    float gWx = 0.0f;
    float gWy = 0.0f;
    float gWz = 0.0f;
    
    float gradrhox = 0.0f;
    float gradrhoy = 0.0f;
    float gradrhoz = 0.0f;
    
    int j = 0;

    for (int jj = 0; jj < MAX_ngb; jj++)
    {
    
      j = ngb[i * MAX_ngb + jj];
    
      if (j != -1) // Note that the Typ[j] == 0 condition was already checked during ngb array construction!
      {
        dx = x[j] - x[i];
        dy = y[j] - y[i];
        dz = z[j] - z[i];

        rr = sqrt(dx * dx + dy * dy + dz * dz);
        hij = 0.5f * (h[i] + h[j]);
        q = rr / hij;

        if (q <= 2.0f)
        {

          nW = 0.0f;
          gWx = 0.0f;
          gWy = 0.0f;
          gWz = 0.0f;

          sig = 1.0f / M_PI;
          hij5 = hij * hij * hij * hij * hij;

          if (q <= 1.0f)
          {
            nW = sig / hij5 * (-3.0f + (9.0f / 4.0f) * q);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          if ((q > 1.0f) && (q <= 2.0f))
          {
            nW = -3.0f * sig / (4.0f * hij5) * (2.0f - q) * (2.0f - q) / (q + 1e-10);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          gradrhox += mass[j] * gWx;
          gradrhoy += mass[j] * gWy;
          gradrhoz += mass[j] * gWz;
        }
      }
    }
    gradRhoNorm[i] = sqrt(gradrhox * gradrhox + gradrhoy * gradrhoy + gradrhoz * gradrhoz);
  }
}


//=====================================================
//================== div_curlVel_ngb ======================
//=====================================================
__global__ void div_curlVel_ngb(int *Typ, float *divV, float *curlV, float *x, float *y, float *z,
                            float *vx, float *vy, float *vz, float *rho, float *mass,
                            int *ngb, int MAX_ngb, float *h, int N)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float dx, dy, dz, rr, hij, q, vxji, vyji, vzji, hij5, sig;
    float nW = 0.0f;
    float gWx = 0.0f;
    float gWy = 0.0f;
    float gWz = 0.0f;
    float vxij, vyij, vzij;
    float ss = 0.0f;
    float curlVx = 0.0f;
    float curlVy = 0.0f;
    float curlVz = 0.0f;
    
    int j = 0;

    for (int jj = 0; jj < MAX_ngb; jj++)
    {
    
      j = ngb[i * MAX_ngb + jj];
    
      if (j != -1) // Note that the Typ[j] == 0 condition was already checked during ngb array construction!
      {
        dx = x[j] - x[i];
        dy = y[j] - y[i];
        dz = z[j] - z[i];

        rr = sqrt(dx * dx + dy * dy + dz * dz);
        hij = 0.5f * (h[i] + h[j]);
        q = rr / hij;

        if (q <= 2.0f)
        {

          nW = 0.0f;
          gWx = 0.0f;
          gWy = 0.0f;
          gWz = 0.0f;

          sig = 1.0f / M_PI;
          hij5 = hij * hij * hij * hij * hij;

          if (q <= 1.0f)
          {
            nW = sig / hij5 * (-3.0f + (9.0f / 4.0f) * q);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          if ((q > 1.0f) && (q <= 2.0f))
          {
            nW = -3.0f * sig / (4.0f * hij5) * (2.0f - q) * (2.0f - q) / (q + 1e-10);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          vxji = vx[j] - vx[i];
          vyji = vy[j] - vy[i];
          vzji = vz[j] - vz[i];

          ss += mass[j] / rho[i] * (vxji * gWx + vyji * gWy + vzji * gWz);

          vxij = vx[i] - vx[j]; //-vxji;
          vyij = vy[i] - vy[j]; //-vyji;
          vzij = vz[i] - vz[j]; //-vzji;

          curlVx += mass[j] / rho[i] * (vyij * gWz - vzij * gWy); // eq. 18 in Beck et al. 2016.
          curlVy += mass[j] / rho[i] * (vzij * gWx - vxij * gWz);
          curlVz += mass[j] / rho[i] * (vxij * gWy - vyij * gWx);
        }
      }
    }
    divV[i] = ss; // abs(ss);
    curlV[i] = sqrt(curlVx * curlVx + curlVy * curlVy + curlVz * curlVz);
  }
}

//===========================================================
//====================== acc_sph_ngb ============================
//===========================================================
__global__ void acc_sph_ngb(int *Typ, float *x, float *y, float *z, float *vx, float *vy, float *vz,
                        float *h, float *c, float *rho, float *divV, float *curlV,
                        float *mass, float *P, float *ax, float *ay, float *az,
                        int *ngb, int MAX_ngb, float visc_alpha, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float dx, dy, dz, rr, hij, q, sig, hij5, gWx, gWy, gWz, nW;
    float vxij, vyij, vzij, wij, vij_rij, vij_sig, rhoij, PIij, fi, fj, fij;
    float axt = 0.0f;
    float ayt = 0.0f;
    float azt = 0.0f;
    
    int j = 0;

    for (int jj = 0; jj < MAX_ngb; jj++)
    {
    
      j = ngb[i * MAX_ngb + jj];
    
      if (j != -1) // Note that the Typ[j] == 0 condition was already checked during ngb array construction!
      {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        dz = z[i] - z[j];

        rr = sqrt(dx * dx + dy * dy + dz * dz);

        hij = 0.5f * (h[i] + h[j]);

        if (rr < 2.0f * hij)
        {

          nW = 0.0f;
          gWx = 0.0f;
          gWy = 0.0f;
          gWz = 0.0f;
          sig = 1.0f / M_PI;
          hij5 = hij * hij * hij * hij * hij;
          q = rr / hij;

          if (q <= 1.0f)
          {
            nW = sig / hij5 * (-3.0f + (9.0f / 4.0f) * q);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          if ((q > 1.0f) && (q <= 2.0f))
          {
            nW = -3.0f * sig / (4.0f * hij5) * (2.0f - q) * (2.0f - q) / (q + 1e-10);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          //-------- PIij ---------
          vxij = vx[i] - vx[j];
          vyij = vy[i] - vy[j];
          vzij = vz[i] - vz[j];

          vij_rij = vxij * dx + vyij * dy + vzij * dz;

          wij = vij_rij / (rr + 1e-5);
          vij_sig = c[i] + c[j] - 3.0f * wij;
          rhoij = 0.5f * (rho[i] + rho[j]);

          PIij = 0.0f;
          if (vij_rij <= 0.0f)
          {

            PIij = -0.5f * visc_alpha * vij_sig * wij / rhoij;

            //------- Shear-viscosity correction -------
            fi = abs(divV[i]) / (abs(divV[i]) + curlV[i] + 0.0001 * c[i] / h[i]);
            fj = abs(divV[j]) / (abs(divV[j]) + curlV[j] + 0.0001 * c[j] / h[j]);
            fij = 0.5f * (fi + fj);
            PIij = fij * PIij;
            //------- End of Shear-visc. correction -----
          }

          axt -= mass[j] * (P[i] / rho[i] / rho[i] + P[j] / rho[j] / rho[j] + PIij) * gWx;
          ayt -= mass[j] * (P[i] / rho[i] / rho[i] + P[j] / rho[j] / rho[j] + PIij) * gWy;
          azt -= mass[j] * (P[i] / rho[i] / rho[i] + P[j] / rho[j] / rho[j] + PIij) * gWz;
        }
      }
    }
    ax[i] = axt;
    ay[i] = ayt;
    az[i] = azt;
  }
}

//===========================================================
//====================== acc_tot ============================
//===========================================================
__global__ void acc_g_sph(int *Typ, float *acc_totx, float *acc_toty, float *acc_totz,
                          float *acc_gx, float *acc_gy, float *acc_gz,
                          float *acc_sphx, float *acc_sphy, float *acc_sphz,
                          int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && ((Typ[i] == 0) || (Typ[i] == 1)))
  {
    acc_totx[i] = acc_gx[i] + acc_sphx[i];
    acc_toty[i] = acc_gy[i] + acc_sphy[i];
    acc_totz[i] = acc_gz[i] + acc_sphz[i];
  }
}

//===============================================
//=================== get_dU_ngb ====================
//===============================================
__global__ void get_dU_ngb(int *Typ, float *x, float *y, float *z, float *vx, float *vy, float *vz,
                       float *h, float *c, float *rho, float *divV, float *curlV,
                       float *mass, float *P, float *dudt,
                       int *ngb, int MAX_ngb, float visc_alpha, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    float dx, dy, dz, rr, hij, q, sig, hij5, gWx, gWy, gWz, nW, vij_gWij;
    float vxij, vyij, vzij, wij, vij_rij, vij_sig, rhoij, PIij, fi, fj, fij;
    float dut = 0.0f;
    
    int j = 0;

    for (int jj = 0; jj < MAX_ngb; jj++)
    {
    
      j = ngb[i * MAX_ngb + jj];
    
      if (j != -1) // Note that the Typ[j] == 0 condition was already checked during ngb array construction!
      {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        dz = z[i] - z[j];

        rr = sqrt(dx * dx + dy * dy + dz * dz);

        hij = 0.5f * (h[i] + h[j]);

        if (rr < 2.0f * hij)
        {

          nW = 0.0f;
          gWx = 0.0f;
          gWy = 0.0f;
          gWz = 0.0f;
          sig = 1.0f / M_PI;
          hij5 = hij * hij * hij * hij * hij;
          q = rr / hij;

          if (q <= 1.0f)
          {
            nW = sig / hij5 * (-3.0f + (9.0f / 4.0f) * q);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          if ((q > 1.0f) && (q <= 2.0f))
          {
            nW = -3.0f * sig / (4.0f * hij5) * (2.0f - q) * (2.0f - q) / (q + 1e-10);
            gWx = nW * dx;
            gWy = nW * dy;
            gWz = nW * dz;
          }

          //-------- PIij ---------
          vxij = vx[i] - vx[j];
          vyij = vy[i] - vy[j];
          vzij = vz[i] - vz[j];

          vij_gWij = vxij * gWx + vyij * gWy + vzij * gWz;

          vij_rij = vxij * dx + vyij * dy + vzij * dz;

          wij = vij_rij / (rr + 1e-5);
          vij_sig = c[i] + c[j] - 3.0f * wij;
          rhoij = 0.5f * (rho[i] + rho[j]);

          PIij = 0.0f;
          if (vij_rij <= 0.0f)
          {

            PIij = -0.5f * visc_alpha * vij_sig * wij / rhoij;

            //------- Shear-viscosity correction -------
            fi = abs(divV[i]) / (abs(divV[i]) + curlV[i] + 0.0001 * c[i] / h[i]);
            fj = abs(divV[j]) / (abs(divV[j]) + curlV[j] + 0.0001 * c[j] / h[j]);
            fij = 0.5f * (fi + fj);
            PIij = fij * PIij;
            //------- End of Shear-visc. correction -----
          }
          dut += mass[j] * (P[i] / rho[i] / rho[i] + PIij / 2.0f) * vij_gWij;
        }
      }
    }
    dudt[i] = dut;
  }
}

//==================================================
//============== update u & utprevious =============
//==================================================
__global__ void u_updater(int *Typ, float *u, float *dudt,
                          float *utprevious, float dt, int N)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    u[i] = u[i] + 0.5f * dt * (dudt[i] + utprevious[i]);
    utprevious[i] = dudt[i];
  }
}

//===========================================================
//================= velocity evolution ======================
//===========================================================
__global__ void v_evolve(int *Typ, float *vx, float *vy, float *vz,
                         float *accx, float *accy, float *accz,
                         float dt, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && ((Typ[i] == 0) || (Typ[i] == 1)))
  {
    vx[i] += accx[i] * dt / 2.0f;
    vy[i] += accy[i] * dt / 2.0f;
    vz[i] += accz[i] * dt / 2.0f;
  }
}

//===========================================================
//================= position evolution ======================
//===========================================================
__global__ void r_evolve(int *Typ, float *x, float *y, float *z,
                         float *vx, float *vy, float *vz,
                         float dt, int ndx_BH, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && ((Typ[i] == 0) || (Typ[i] == 1)))
  {
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
    z[i] += vz[i] * dt;
    
    if (i == ndx_BH)
    {
      // Re-setting the BH position back to (0, 0, 0)!
      x[ndx_BH] = 0.0f;
      y[ndx_BH] = 0.0f;
      z[ndx_BH] = 0.0f;
    }
  }
}


//===========================================================
//================= position evolution ======================
//===========================================================
__global__ void r_evolveT(int *Typ, float *x, float *y, float *z,
                         float *vx, float *vy, float *vz,
                         float dt, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && ((Typ[i] == 0) || (Typ[i] == 1)))
  {
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
    z[i] += vz[i] * dt;
  }
}

//===========================================================
//=================== dt estimation =========================
//===========================================================
__global__ void dt_array_indiv_dt(int *Typ, float *x, float *y, float *z,
                                  float *vx, float *vy, float *vz,
                                  float *accx, float *accy, float *accz,
                                  float *accx_tot, float *accy_tot, float *accz_tot,
                                  float *h, float *c, float *dt_particles,
                                  float *abs_acc_g, float *abs_acc_tot,
                                  float *divV, float *dh_dt, float C_CFL,
                                  float visc_alpha, float *eps, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {

    abs_acc_g[i] = sqrt(accx[i] * accx[i] + accy[i] * accy[i] + accz[i] * accz[i]);
    abs_acc_tot[i] = sqrt(accx_tot[i] * accx_tot[i] + accy_tot[i] * accy_tot[i] + accz_tot[i] * accz_tot[i]);

    float dx, dy, dz, vxij, vyij, vzij, wij, vij_rij, vij_sig, rr;

    float max_vij_sig = 0.0f;
    float tmp = 0.0f;
    for (int j = 0; j < N; j++)
    {
      if (Typ[j] == 0)
      {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        dz = z[i] - z[j];

        rr = sqrt(dx * dx + dy * dy + dz * dz);

        vxij = vx[i] - vx[j];
        vyij = vy[i] - vy[j];
        vzij = vz[i] - vz[j];

        vij_rij = vxij * dx + vyij * dy + vzij * dz;

        wij = vij_rij / (rr + 1e-5);

        vij_sig = c[i] + c[j] - 3.0f * wij;

        tmp = vij_sig;

        if (tmp > max_vij_sig)
        {
          max_vij_sig = tmp;
        }
      }
    }

    float dt_hyd = C_CFL * h[i] / max_vij_sig; // eq. 16 in Springel et al - 2005.

    // dt_cour: ref: Gadget 2 paper, eq. 49.
    float abs_v_i = sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    float dt_cour = C_CFL * h[i] / (h[i] * abs(divV[i]) + max(c[i], abs_v_i) * (1.0f + 0.6f * visc_alpha));

    // dt_G1, dt_G2 ref: Dave et al - 1997 (Parallel TreeSPH) sec. 3.7.
    float ettaa = 0.4;
    float dt_G1 = ettaa * sqrt(eps[i] / abs_acc_tot[i]);
    float dt_G2 = ettaa * eps[i] / abs_v_i;

    float dt_f = sqrt(h[i] / abs_acc_g[i]);
    float dt_kin = sqrt(h[i] / abs_acc_tot[i]);

    float dtZ[6] = {dt_hyd, dt_cour, dt_f, dt_kin, dt_G1, dt_G2}; // Note to modify the max range in for loop
    float dtxx = dtZ[0];
    // finding the minimum !!
    for (int k = 0; k < 6; k++) // modify if you changed dtZ[]!!! IMPORTANT!!
    {
      if (dtZ[k] <= dtxx)
      {
        dtxx = dtZ[k];
      }
    }
    dt_particles[i] = dtxx;

    dh_dt[i] = 1.0f / 3.0f * h[i] * divV[i]; // See the line below eq.31 in Gadget 2 paper.
  }
}



//==============================================================
//====================== hash function =========================
//==============================================================
__device__ unsigned long long hash(unsigned long long x, unsigned long long y) {
    return x * 6364136223846793005ULL + y;
}



//==============================================================
//========= Outflow injection 2 (Richings et al - 2018)===========
//==============================================================
__global__ void outflow_injector2(int *Typ, float *x, float *y, float *z,
                                 float *vx, float *vy, float *vz,
                                 float *h, float *eps, float *mass,
                                 float Nngb_f, float *Nngb_previous,
                                 float *u, float M_dot_in, float v_in,
                                 float mass_sph_high_res, float u_for_10K_Temp,
                                 float *leftover_mass, float dt, int ndx_BH, int N,
                                 unsigned long long seed)
{

  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if (idx == 0)
  {
    // Calculate the mass to inject in each time step
    float mass_to_inject = M_dot_in * dt;
    mass_to_inject += *leftover_mass;

    // Calculate number of particle pairs to inject
    int num_pairs = static_cast<int>(mass_to_inject / (2.0f * mass_sph_high_res));

    float one_pc_in_code_unit = 0.001f; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    float h_outflow = 0.001; // ONLY place holder!

    // Finding the last index for which Typ is not equal to -1.
    // Note that with the way I do in this for loop, jj will be the number of times Typ is not
    // equal to -1. So to convert it to index, we need to subtract 1 from it !!!
    int jj = 0;
    for (int i = 0; i < N; i++)
    {
      if ((Typ[i] == 0) || (Typ[i] == 1))
      {
        jj += 1;
      }
    }
    jj -= 1; // See the above explanation!!
    
    //--- Estimation of h_outflow --
    if (jj > ndx_BH)
    {
      h_outflow = h[jj]; // We take the h of the last injected outflow particle!
    }
    else
    {
      h_outflow = h[ndx_BH - 1]; // This happens when we had not yet injected any outflow particle!
    }

    // Spawn particle pairs
    for (int j = 0; j < num_pairs; j++)
    {

      // Initialize cuRAND generator with the provided seed
      curandState state;
      curand_init(hash(seed, j), idx, 0, &state);

      // Calculate random distance from 0 to 1 pc (in cm)
      float rt = curand_uniform(&state) * one_pc_in_code_unit;
      // Calculate random orientation (assuming spherical coordinates)
      float theta = curand_uniform(&state) * 2.0f * M_PI;
      //float phi = curand_uniform(&state) * M_PI; // This was non symmetric!!
      float phi = acos(2.0f * curand_uniform(&state) - 1.0f);
      // Calculate Cartesian coordinates of the pair
      float xt = rt * sin(phi) * cos(theta);
      float yt = rt * sin(phi) * sin(theta);
      float zt = rt * cos(phi);

      float rr = sqrt(xt * xt + yt * yt + zt * zt);

      float vxt = xt / rr * v_in;
      float vyt = yt / rr * v_in;
      float vzt = zt / rr * v_in;

      //---- Injecting the first particle of the pair ----
      Typ[jj + 1] = 0;

      x[jj + 1] = xt;
      y[jj + 1] = yt;
      z[jj + 1] = zt;

      vx[jj + 1] = vxt;
      vy[jj + 1] = vyt;
      vz[jj + 1] = vzt;

      mass[jj + 1] = mass_sph_high_res;
      h[jj + 1] = h_outflow;
      eps[jj + 1] = h_outflow;
      u[jj + 1] = u_for_10K_Temp;

      Nngb_previous[jj + 1] = Nngb_f;

      //---- Injecting the second particle of the pair ----
      Typ[jj + 2] = 0;

      x[jj + 2] = -xt;
      y[jj + 2] = -yt;
      z[jj + 2] = -zt;

      vx[jj + 2] = -vxt;
      vy[jj + 2] = -vyt;
      vz[jj + 2] = -vzt;

      mass[jj + 2] = mass_sph_high_res;
      h[jj + 2] = h_outflow;
      eps[jj + 2] = h_outflow;
      u[jj + 2] = u_for_10K_Temp;

      Nngb_previous[jj + 2] = Nngb_f;

      jj += 2;

      mass_to_inject -= 2.0f * mass_sph_high_res;      
    }
    *leftover_mass = mass_to_inject;
  }
}


//==================================================================================
//========= Isothermal gravitational acceleration (Richings et al - 2018)===========
//==================================================================================
__global__ void galaxy_isothermal_potential(int *Typ, float *x, float *y, float *z, float *accx,
                                            float *accy, float *accz, float sigma, float G, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && ((Typ[i] == 0) || (Typ[i] == 1)))
  {

    float rr = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

    // Compute the enclosed mass at radius rr (eq.1 in Richings et al - 2018).
    float enclosed_mass = 2.0 * sigma * sigma * rr / G;

    // Compute the gravitational acceleration
    accx[i] += -G * enclosed_mass * x[i] / rr / rr / rr; // ===> -G * M * x / r^3
    accy[i] += -G * enclosed_mass * y[i] / rr / rr / rr;
    accz[i] += -G * enclosed_mass * z[i] / rr / rr / rr;
  }
}



//=====================================================================================
//======================= Heating, Cooling using CLOUDY =============================== (24 Feb 2024)
//=====================================================================================
//===== findClosestIndex (Is Used for Heating, Cooling)
__device__ int find_nx0(float *A, int N, float x)
{
  int nxClosest = 0; // Start with the first element as the closest
  float smallestDiff = fabsf(A[0] - x); // Initialize the smallest difference using fabsf for float absolute value in CUDA

  for (int i = 1; i < N; ++i)
  {
    float currentDiff = fabsf(A[i] - x); // Compute the current difference using fabsf
    if (currentDiff < smallestDiff) { // If the current difference is smaller, update the closest element
      smallestDiff = currentDiff;
      nxClosest = i;
    }
  }

  int nx0 = 0; // just a place holder.

  if (x >= A[nxClosest])
    nx0 = nxClosest;
  else
    nx0 = nxClosest - 1;

  if (nx0 == N - 1) // To handle cases in which the index is at the last element
    nx0 = nxClosest - 1;

  if (nx0 == -1) // Handle case where computed index is below 0
    nx0 = 0;

  return nx0; // Return nx0. nx1 will be computed elsewhere!
}


//===== interpolate_4d_hypercube (Is Used for Heating, Cooling)
__device__ float interpolate_4d_hypercube(float nH_p, float T_p, float r_p, float NH_p, float *nH, float *T, float *r, float *NH,
                                          float *Gam, float *Lam, int nxnH0, int nxT0, int nxr0, int nxNH0, int nxnH1, int nxT1, int nxr1, int nxNH1,
                                          int N_nH, int N_T, int N_r, int N_NH)
{
  float dx = (nH_p - nH[nxnH0]) / (nH[nxnH1] - nH[nxnH0]);
  float dy = (T_p - T[nxT0]) / (T[nxT1] - T[nxT0]);
  float dz = (r_p - r[nxr0]) / (r[nxr1] - r[nxr0]);
  float dw = (NH_p - NH[nxNH0]) / (NH[nxNH1] - NH[nxNH0]);

  float pG[16];
  float pL[16];

  int stride_nH = N_T * N_r * N_NH;
  int stride_T = N_r * N_NH;
  int stride_r = N_NH;

  int nxs = 0;
  for (int i = nxnH0; i <= nxnH1; i++)
  {
    for (int j = nxT0; j <= nxT1; j++)
    {
      for (int k = nxr0; k <= nxr1; k++)
      {
        for (int l = nxNH0; l<= nxNH1; l++)
        {
          pG[nxs] = Gam[i * stride_nH + j * stride_T + k * stride_r + l]; // heating!
          pL[nxs] = Lam[i * stride_nH + j * stride_T + k * stride_r + l]; // cooling!
          nxs += 1;
        }
      }
    }
  }

  //-------------> HEATING <---------------
  // Interpolate along the x-axis
  float c00 = pG[0] * (1.0f - dx) + pG[8] * dx;
  float c01 = pG[1] * (1.0f - dx) + pG[9] * dx;
  float c10 = pG[2] * (1.0f - dx) + pG[10] * dx;
  float c11 = pG[3] * (1.0f - dx) + pG[11] * dx;
  float c20 = pG[4] * (1.0f - dx) + pG[12] * dx;
  float c21 = pG[5] * (1.0f - dx) + pG[13] * dx;
  float c30 = pG[6] * (1.0f - dx) + pG[14] * dx;
  float c31 = pG[7] * (1.0f - dx) + pG[15] * dx;

  // Correct interpolation along the y-axis
  float c0 = c00 * (1.0f - dy) + c10 * dy;
  float c1 = c01 * (1.0f - dy) + c11 * dy;
  float c2 = c20 * (1.0f - dy) + c30 * dy;
  float c3 = c21 * (1.0f - dy) + c31 * dy;

  // Interpolate along the z-axis
  float c = c0 * (1.0f - dz) + c2 * dz;
  float d = c1 * (1.0f - dz) + c3 * dz;

  // Finally, interpolate along the w-axis
  float GamH = c * (1.0f - dw) + d * dw;

  //-------------> COOLING <---------------
  // Interpolate along the x-axis
  c00 = pL[0] * (1.0f - dx) + pL[8] * dx;
  c01 = pL[1] * (1.0f - dx) + pL[9] * dx;
  c10 = pL[2] * (1.0f - dx) + pL[10] * dx;
  c11 = pL[3] * (1.0f - dx) + pL[11] * dx;
  c20 = pL[4] * (1.0f - dx) + pL[12] * dx;
  c21 = pL[5] * (1.0f - dx) + pL[13] * dx;
  c30 = pL[6] * (1.0f - dx) + pL[14] * dx;
  c31 = pL[7] * (1.0f - dx) + pL[15] * dx;

  // Correct interpolation along the y-axis
  c0 = c00 * (1.0f - dy) + c10 * dy;
  c1 = c01 * (1.0f - dy) + c11 * dy;
  c2 = c20 * (1.0f - dy) + c30 * dy;
  c3 = c21 * (1.0f - dy) + c31 * dy;

  // Interpolate along the z-axis
  c = c0 * (1.0f - dz) + c2 * dz;
  d = c1 * (1.0f - dz) + c3 * dz;

  // Finally, interpolate along the w-axis
  float LamC = c * (1.0f - dw) + d * dw;

  return powf(10, GamH) - powf(10, LamC); // Use powf for floating-point exponentiation
}


//===== getGamMinusLam (Is Used for Heating, Cooling)
__device__ float getGamMinusLam(float nH_p, float u_p, float r_p, float NH_p, float *nH, float *T, float *r, float *NH, float *uarr, float *muarr,
                                float *Gam, float *Lam, int N_nH, int N_T, int N_r, int N_NH, float gamma, float kB, float mH,
                                float nHLowBound, float nHUpBound, float TLowBound, float TUpBound, float rLowBound, float rUpBound,
                                float NHLowBound, float NHUpBound)
{
  // Adjusting input parameters to be within the bounds
  nH_p = max(nHLowBound + 0.0001f, min(nH_p, nHUpBound - 0.0001f));
  r_p = max(rLowBound + 0.0001f, min(r_p, rUpBound - 0.0001f));
  NH_p = max(NHLowBound + 0.0001f, min(NH_p, NHUpBound - 0.0001f));

  // Finding indexes for interpolation
  int nxnH0 = find_nx0(nH, N_nH, nH_p);
  int nxr0 = find_nx0(r, N_r, r_p);
  int nxNH0 = find_nx0(NH, N_NH, NH_p);

  // Convert u_p to T_p
  u_p = log10(u_p);

  int S1 = N_T * N_r * N_NH;
  int S2 = N_r * N_NH;
  int S3 = N_NH;

  int nxnH1 = nxnH0 + 1;
  int nxr1 = nxr0 + 1;
  int nxNH1 = nxNH0 + 1;

  int ndx = 0;
  int nxu0 = -1;
  for (int j = 0; j < N_T; j++)
  {
    ndx = (nxnH1 * S1) + (j * S2) + (nxr1 * S3) + nxNH1;
    if (nxu0 == -1 && u_p <= uarr[ndx])
    {
      nxu0 = ndx;
      break; // Found the first index where u_p <= uarr[ndx], stop searching
    }
  }

  // Handling the case where u_p is greater than max(uarr)
  if (nxu0 == -1)
    nxu0 = N_nH * N_T * N_r * N_NH - 1;

  float mu_p = muarr[nxu0];
  float T_p = (gamma - 1.0f) * mu_p * mH / kB * powf(10, u_p);
  T_p = log10(T_p);

  // Ensuring T_p is within bounds
  T_p = max(TLowBound + 0.0001f, min(T_p, TUpBound - 0.0001f));

  int nxT0 = find_nx0(T, N_T, T_p);
  int nxT1 = nxT0 + 1;

  // Performing 4D interpolation
  float Gam_minus_Lam = interpolate_4d_hypercube(nH_p, T_p, r_p, NH_p, nH, T, r, NH, Gam, Lam, nxnH0, nxT0, nxr0, nxNH0,
                                                 nxnH1, nxT1, nxr1, nxNH1, N_nH, N_T, N_r, N_NH);

  return Gam_minus_Lam;
}


//===== doCooling (Is Used for Heating, Cooling)
__global__ void doCooling(int *Typ, float *x, float *y, float *z, float *rho_p, float *u_p, float *NHtot,
                          float *nHGrid, float *TGrid, float *rGrid, float *NHGrid, float *uarr, float *muarr, float *Gam, float *Lam, 
                          int N_nH, int N_T, int N_r, int N_NH, float gamma, float kB, float mH,
                          float nHLowBound, float nHUpBound, float TLowBound, float TUpBound, float rLowBound, float rUpBound,
                          float NHLowBound, float NHUpBound, float dt,
                          float unit_time_in_sec, float unit_density_in_cgs, float unit_length_in_kpc, float unit_u_in_cgs, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  
  if ((i < N) && (Typ[i] == 0))
  {
    float XH = 0.7;

    float dt_s = dt * unit_time_in_sec;

    float rho_cgs = rho_p[i] * unit_density_in_cgs;
    float nH_p = XH * rho_cgs / mH;
    nH_p = log10(nH_p);

    float u_ad = u_p[i] * unit_u_in_cgs; // log10 of it will be taken inside the "getGamMinusLam" function!

    float r_p = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
    r_p = r_p * unit_length_in_kpc;

    float NH_p = NHtot[i]; // !!!!!! MAKE sure that NHtot are already in log10 form!!!!!!!!!!!!

    float GL = 0.0;
    

    float u_old = u_ad; 

    float u = u_old;
    float u_upper = u;

    // log10 of u will be taken inside the "getGamMinusLam" function!
    GL = getGamMinusLam(nH_p, u, r_p, NH_p, nHGrid, TGrid, rGrid, NHGrid, uarr, muarr,
                        Gam, Lam, N_nH, N_T, N_r, N_NH, gamma, kB, mH,
                        nHLowBound, nHUpBound, TLowBound, TUpBound, rLowBound, rUpBound, NHLowBound, NHUpBound);
    float f_u = u - (u_ad + (GL/rho_cgs*dt_s));

    float u_lower = u; // just a place holder!
    float f_u_upper = f_u;

    //-----> (f_u > 0.0) <--------
    if (f_u > 0.0)
    {
      while (f_u_upper > 0.0)
      {
        GL = getGamMinusLam(nH_p, u_upper, r_p, NH_p, nHGrid, TGrid, rGrid, NHGrid, uarr, muarr,
                            Gam, Lam, N_nH, N_T, N_r, N_NH, gamma, kB, mH,
                            nHLowBound, nHUpBound, TLowBound, TUpBound, rLowBound, rUpBound, NHLowBound, NHUpBound);
        f_u_upper = u_upper - (u_ad + (GL/rho_cgs*dt_s));
        u_upper = u_upper / 1.1;
      }
        
      u_lower = u_upper;
      u_upper = u_lower * 1.1;
    }

    //-----> (f_u < 0.0) <--------
    if (f_u < 0.0)
    {    
      while (f_u_upper < 0.0)
      {
        GL = getGamMinusLam(nH_p, u_upper, r_p, NH_p, nHGrid, TGrid, rGrid, NHGrid, uarr, muarr,
                            Gam, Lam, N_nH, N_T, N_r, N_NH, gamma, kB, mH,
                            nHLowBound, nHUpBound, TLowBound, TUpBound, rLowBound, rUpBound, NHLowBound, NHUpBound);
        f_u_upper = u_upper - (u_ad + (GL/rho_cgs*dt_s));
        u_upper = u_upper * 1.1;
      }
      u_upper = u_upper;
      u_lower = u_upper / 1.1;
    }

    int MAXniter = 100;
    int niter = 0;
    float du = u;

    while (du/u > 1e-4  && niter < MAXniter)
    {
      u = 0.5f * (u_lower + u_upper);
      
      GL = getGamMinusLam(nH_p, u, r_p, NH_p, nHGrid, TGrid, rGrid, NHGrid, uarr, muarr,
                          Gam, Lam, N_nH, N_T, N_r, N_NH, gamma, kB, mH,
                          nHLowBound, nHUpBound, TLowBound, TUpBound, rLowBound, rUpBound, NHLowBound, NHUpBound);
      f_u = u - (u_ad + (GL/rho_cgs*dt_s));
      
      if (f_u > 0.0)
        u_upper = u;
      else
        u_lower = u;
      
      du = abs(u_upper - u_lower);
      
      niter++;
      
      if (niter == 100)
        printf("MAXniter in doCooling function REACHED !!!!!!!!!!\n");
    }
    
    u_p[i] = u / unit_u_in_cgs;
    
  }
}



//=====================================================================================
//========================== For debugging Cooling ====================================
//=====================================================================================
__global__ void u_before_adiabatic(int *Typ, float *u, float *uB, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    uB[i] = u[i];
  }

}


//=====================================================================================
//========================== For debugging Cooling ====================================
//=====================================================================================
__global__ void u_After_cooling(int *Typ, float *u, float *uB, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    uB[i] = u[i];
  }

}



//=====================================================================================
//========================== For debugging Cooling ====================================
//=====================================================================================
__global__ void dudt_previous(int *Typ, float *u, float *uB, int N)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if ((i < N) && (Typ[i] == 0))
  {
    uB[i] = u[i];
  }

}


#endif
