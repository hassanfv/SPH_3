
# In this version, I include grain_recombination (29 May 2024).

#import h5py
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from data5 import *
from scipy.interpolate import RegularGridInterpolator


# Cooling due to free-free emission ::: will be multiplied by (nHp+nHep+nHepp)*ne later in the code.
def gfree(T):

  gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2/3.0)
  g4_val = 1.42e-27 * gff * T**0.5
  
  return g4_val



#----- Lambda
def Lambda(T, nH0, nHp, nHe0, nHep, nHepp, nC0, nC1, nC2, nC3, nC4, nC5, nC6):
    Tx = np.log10(T)
    ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
    
    cFree = nHp + nHep + 4.0 * nHepp + nC1 + 4.0*nC2 + 9.*nC3 + 16. * nC4 + 25.0 * nC5 + 36.0 * nC6
    
    #----- # Glover & Jappsen - 2007 -----
    z = 0.0 # current time redshift!
    TCMB_0 = 2.7255
    TCMB = TCMB_0 * (1.0 + z)
    LCompton = 1.017e-37 * TCMB**4 * (T - TCMB) * ne
    #--------------------------------------
    
    #------> Heating by Cosmic Rays (see eq_9 in Goldsmith & Langer - 1978) <--------
    dQcr = 20. * 1.60218e-12 # erg
    
    Lcr_H0  = 1.00 * cr_rateH * dQcr * nH0  # erg/cm^-3/s
    Lcr_He0 = 1.09 * cr_rateH * dQcr * nHe0
    Lcr_Hep = 0.25 * cr_rateH * dQcr * nHep
    Lcr_C0  = 3.83  * cr_rateH * dQcr * nC0 # erg/cm^-3/s
    Lcr_C1  = 1.66  * cr_rateH * dQcr * nC1
    Lcr_C2  = 0.83  * cr_rateH * dQcr * nC2
    Lcr_C3  = 0.45  * cr_rateH * dQcr * nC3
    Lcr_C4  = 0.069 * cr_rateH * dQcr * nC4
    Lcr_C5  = 0.028 * cr_rateH * dQcr * nC5
    
    grain_cool_H_1_0 = grain_cool_rate(T, ne, nH0, nHp, G0, A_v, dust_ratio, Temp, Psi)
    grain_cool_He_1_0 = grain_cool_rate(T, ne, nH0, nHp, G0, A_v, dust_ratio, Temp, Psi)
    grain_cool_C_1_0 = grain_cool_rate(T, ne, nH0, nHp, G0, A_v, dust_ratio, Temp, Psi)
    
    Lamb = (
          10**g1(Tx) * ne * nH0  # H0
        + 10**g2(Tx) * ne * nHp # Hp
        + 10**grain_cool_H_1_0 * nHp * ne
        + 10**g3(Tx) * nHe0 * ne # He0 
        + 10**g4(Tx) * nHep * ne # Hep
        + 10**grain_cool_He_1_0 * nHep * ne
        + 10**g5(Tx) * nHepp * ne# Hepp
        + 10**C0_cooling_rate(T, nH0, ne, nHp, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nC0 * ne # cooling via C0
        + 10**Cp_cooling_rate(T, ne, Temp_2d, elecDensity_2d) * nC1 * ne # cooling via Cp or C1
        + 10**grain_cool_C_1_0 * nC1 * ne
        + 10**gC2(Tx) * nC2 * ne # C2
        + 10**gC3(Tx) * nC3 * ne # C3
        + 10**gC4(Tx) * nC4 * ne # C4
        + 10**gC5(Tx) * nC5 * ne # C5
        + 10**gC6(Tx) * nC6 * ne # C6
        + gfree(T) * ne * cFree # free-free emission
        + LCompton
        - (Lcr_H0 + Lcr_He0 + Lcr_Hep + Lcr_C0 + Lcr_C1 + Lcr_C2 + Lcr_C3 + Lcr_C4 + Lcr_C5) # minus will be multiplied with a minus and become + !!!!
    )
    return Lamb



#----- func
def func(t, y):

    nH0, nHp, nHe0, nHep, nHepp, nC0, nC1, nC2, nC3, nC4, nC5, nC6, T = y
    Tx = np.log10(T)
    
    ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
    ntot = ne + (nH0 + nHp) + (nHe0 + nHep + nHepp) + (nC0 + nC1 + nC2 + nC3 + nC4 + nC5 + nC6)
    
    #------> Cosmic Ray rates (see eq_9 in Goldsmith & Langer - 1978) <--------
    cr_rateH = 1.8e-16 # s^-1 # from chimes param file ----> grid_noneq_evolution.param!
    
    Rcr_H0  = 1.00 * cr_rateH
    Rcr_He0 = 1.09 * cr_rateH
    Rcr_Hep = 0.25 * cr_rateH
    Rcr_C0  = 3.83  * cr_rateH
    Rcr_C1  = 1.66  * cr_rateH
    Rcr_C2  = 0.83  * cr_rateH
    Rcr_C3  = 0.45  * cr_rateH
    Rcr_C4  = 0.069 * cr_rateH
    Rcr_C5  = 0.028 * cr_rateH
    
    grain_rec_H_1_0 =  grain_recomb_rate('Hp', T, ne, G0, A_v, Temp, Psi)
    grain_rec_He_1_0 = grain_recomb_rate('Hep', T, ne, G0, A_v, Temp, Psi)
    grain_rec_C_1_0 =  grain_recomb_rate('C1', T, ne, G0, A_v, Temp, Psi)
    
    dnH0_dt = 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne - Rcr_H0 * nH0 + 10**grain_rec_H_1_0 * nHp * ne
    dnHp_dt = 10**k1(Tx) * nH0 * ne - 10**k2(Tx) * nHp * ne + Rcr_H0 * nH0 - 10**grain_rec_H_1_0 * nHp * ne
    
    dnHe0_dt = 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne - Rcr_He0 * nHe0 + 10**grain_rec_He_1_0 * nHep * ne
    
    dnHep_dt = (10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne +
               Rcr_He0 * nHe0 - Rcr_Hep * nHep - 10** grain_rec_He_1_0 * nHep * ne)

    dnHepp_dt = 10**k4(Tx) * nHep * ne - 10**k6(Tx) * nHepp * ne + Rcr_Hep * nHep
    
    dnC0_dt = (10**C_1_0(Tx) * ne * nC1 - 10**C_0_1(Tx) * ne * nC0
             - 10**ct_C_0_1_Hep(Tx) * nC0 * nHep - 10**ct_C_0_1_Hp(Tx) * nC0 * nHp
             + 10**ct_C_1_0_H0(Tx) * nC1 * nH0 - Rcr_C0 * nC0 + 10**grain_rec_C_1_0 * nC1 * ne) 
    
    dnC1_dt = (10**C_0_1(Tx) * ne * nC0 + 10**C_2_1(Tx) * ne * nC2 - 10**C_1_0(Tx) * ne * nC1 - 10**C_1_2(Tx) * ne * nC1
             + 10**ct_C_0_1_Hep(Tx) * nC0 * nHep + 10**ct_C_0_1_Hp(Tx) * nC0 * nHp - 10**ct_C_1_0_H0(Tx) * nC1 * nH0
             - 10**ct_C_1_2_Hep(Tx) * nC1 * nHep + 10**ct_C_2_1_H0(Tx) * nC2 * nH0 + Rcr_C0 * nC0 - Rcr_C1 * nC1
             - 10**grain_rec_C_1_0 * nC1 * ne)
    
    dnC2_dt = (10**C_1_2(Tx) * ne * nC1 + 10**C_3_2(Tx) * ne * nC3 - 10**C_2_1(Tx) * ne * nC2 - 10**C_2_3(Tx) * ne * nC2
             + 10**ct_C_1_2_Hep(Tx) * nC1 * nHep - 10**ct_C_2_1_H0(Tx) * nC2 * nH0 + 10**ct_C_3_2_He0(Tx) * nC3 * nHe0
             + 10**ct_C_3_2_H0(Tx) * nC3 * nH0 + Rcr_C1 * nC1 - Rcr_C2 * nC2)
    
    dnC3_dt = (10**C_2_3(Tx) * ne * nC2 + 10**C_4_3(Tx) * ne * nC4 - 10**C_3_2(Tx) * ne * nC3 - 10**C_3_4(Tx) * ne * nC3
             - 10**ct_C_3_2_He0(Tx) * nC3 * nHe0 - 10**ct_C_3_2_H0(Tx) * nC3 * nH0 + 10**ct_C_4_3_H0(Tx) * nC4 * nH0
             + 10**ct_C_4_3_He0(Tx) * nC4 * nHe0 + Rcr_C2 * nC2 - Rcr_C3 * nC3)
    
    
    dnC4_dt = (10**C_3_4(Tx) * ne * nC3 + 10**C_5_4(Tx) * ne * nC5 - 10**C_4_3(Tx) * ne * nC4 - 10**C_4_5(Tx) * ne * nC4
             - 10**ct_C_4_3_H0(Tx) * nC4 * nH0 - 10**ct_C_4_3_He0(Tx) * nC4 * nHe0 + 10**ct_C_5_4_H0(Tx) * nC5 * nH0
             + Rcr_C3 * nC3 - Rcr_C4 * nC4)
    
    dnC5_dt = (10**C_4_5(Tx) * ne * nC4 + 10**C_6_5(Tx) * ne * nC6 - 10**C_5_4(Tx) * ne * nC5 - 10**C_5_6(Tx) * ne * nC5
             - 10**ct_C_5_4_H0(Tx) * nC5 * nH0 + Rcr_C4 * nC4 - Rcr_C5 * nC5)
    
    dnC6_dt = 10**C_5_6(Tx) * ne * nC5 - 10**C_6_5(Tx) * ne * nC6 + Rcr_C5 * nC5
    
    Lamb = Lambda(T, nH0, nHp, nHe0, nHep, nHepp, nC0, nC1, nC2, nC3, nC4, nC5, nC6)
    
    dne_dt = dnHp_dt + (dnHep_dt + 2.0 * dnHepp_dt) + (dnC1_dt + 2.0 * dnC2_dt + 3.0 * dnC3_dt + 4.0 * dnC4_dt + 5.0 * dnC5_dt + 6.0 * dnC6_dt)
    dntot_dt = dne_dt + (dnH0_dt + dnHp_dt) + (dnHe0_dt + dnHep_dt + dnHepp_dt) + (dnC0_dt + dnC1_dt + dnC2_dt + dnC3_dt + dnC4_dt + dnC5_dt + dnC6_dt)
    
    dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * (Lamb + 1. / (gamma - 1.) * kB * T * dntot_dt)

    return [dnH0_dt, dnHp_dt, dnHe0_dt, dnHep_dt, dnHepp_dt, dnC0_dt, dnC1_dt, dnC2_dt, dnC3_dt, dnC4_dt, dnC5_dt, dnC6_dt, dT_dt]


cr_rateH = 1.8e-16 # s^-1 # from chimes param file ----> grid_noneq_evolution.param!


nH = 1000.0

He_solar = 10**(-1.0) # See Table_1 in Wiersma et al - 2009, 393, 99–107
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

C_solar = 10**(-3.61) # See Table_1 in Wiersma et al - 2009, 393, 99–107
nC = C_solar * nH

print('nC (cm^-3) = ', nC)

print()
print('H, C before = ', nH, nC)

T_i = 10**6.00

nH0_i = 0.001 * nH
nHp_i = nH - nH0_i

nHe0_i = 0.0001 * nHe
nHep_i = 0.001 * nHe
nHepp_i= nHe - nHe0_i - nHep_i
 
nC0_i = 1e-5 * nC
nC1_i = 1e-5 * nC
nC2_i = 1e-4 * nC
nC3_i = 1e-4 * nC
nC4_i = 1e-3 * nC
nC5_i = 1e-2 * nC
nC6_i = nC - (nC0_i + nC1_i + nC2_i + nC3_i + nC4_i + nC5_i)

y0 = [nH0_i, nHp_i, nHe0_i, nHep_i, nHepp_i, nC0_i, nC1_i, nC2_i, nC3_i, nC4_i, nC5_i, nC6_i, T_i]


A_v = 1.0
G0 = 0.01
dust_ratio = 1.0

t_span = (1*3.16e7, 3000*3.16e7)

#solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True, max_step = 3.16e6)
solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)
#solution = solve_ivp(func, t_span, y0, method='LSODA', jac = jacobian, dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7



nH0  = y[0, :]
nHp  = y[1, :]
nHe0 = y[2, :]
nHep = y[3, :]
nHepp= y[4, :]

nC0 = y[5, :]
nC1 = y[6, :]
nC2 = y[7, :]
nC3 = y[8, :]
nC4 = y[9, :]
nC5 = y[10, :]
nC6 = y[11, :] 
T = y[12, :]

print('nC0 = ', nC0)
print()
print('nC1 = ', nC1)
print()
print('nC2 = ', nC2)
print()
print('nC3 = ', nC3)
print()
print('nC4 = ', nC4)
print()
print('nC5 = ', nC5)
print()
print('nC6 = ', nC6)
print()



#------ Result from "test_primordial_hdf5_v2.py" code -----
with open('chimesRes_C.pkl', 'rb') as f:
  df = pickle.load(f)
# dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp}
t_Arr_in_yrsx = df['t_Arr_in_yrs']
TEvolx = df['TEvol']
nH0x = df['nH0']
nHpx = df['nHp']
nHx = nH0x + nHpx
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx

nC0x = df['nC0']
nC1x = df['nC1']
nC2x = df['nC2']
nC3x = df['nC3']
nC4x = df['nC4']
nC5x = df['nC5']
nC6x = df['nC6']
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x
#----------------------------------------------------------



ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
nex = nHpx + (nHepx + 2.0 * nHeppx) + (nC1x + 2.0 * nC2x + 3.0 * nC3x + 4.0 * nC4x + 5.0 * nC5x + 6.0 * nC6x)



plt.figure(figsize = (16, 8))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs, np.log10(T), s = 2, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 3000)
plt.ylim(1, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.plot(t_yrs, nHe0, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp, color = 'b', label = 'nHepp')

plt.plot(t_Arr_in_yrsx, nHe0x, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx, color = 'b', label = 'nHepp - chimes', linestyle = ':')

plt.xlim(0, 5000)
plt.ylim(1e-8, 300)

plt.yscale('log')
plt.title('solve_ivp')
plt.legend()


plt.subplot(2, 3, 3)
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nH0/nH, label = 'nH0', color = 'r')
plt.plot(T, nHp/nH, label = 'nHp', color = 'g')
plt.plot(T, nHe0/nHeTot, label = 'nHe0', color = 'b')
plt.plot(T, nHep/nHeTot, label = 'nHep', color = 'orange')
plt.plot(T, nHepp/nHeTot,label = 'nHepp', color = 'purple')

plt.plot(TEvolx, nH0x/nHx, label = 'nH0 - chimes', color = 'r', linestyle = ':')
plt.plot(TEvolx, nHpx/nHx, label = 'nHp - chimes', color = 'g', linestyle = ':')
plt.plot(TEvolx, nHe0x/nHeTotx, label = 'nHe0 - chimes', color = 'b', linestyle = ':')
plt.plot(TEvolx, nHepx/nHeTotx, label = 'nHep - chimes', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nHeppx/nHeTotx,label = 'nHepp - chimes', color = 'purple', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)
plt.legend()


ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
#nCC = (nC0 + nC1 + nC2 + nC3 + nC4 + nC5 + nC6)
print('nH after = ', nH)
print('ne after (one ne for each T) = ', ne)
#print('nC after (one nC for each T) = ', nCC)

plt.subplot(2, 3, 4)
plt.plot(T, nC0/nC, label = 'nC0', color = 'r')
plt.plot(T, nC1/nC, label = 'nC1', color = 'g')
plt.plot(T, nC2/nC, label = 'nC2', color = 'b')
plt.plot(T, nC3/nC, label = 'nC3', color = 'orange')
plt.plot(T, nC4/nC, label = 'nC4', color = 'purple')
plt.plot(T, nC5/nC, label = 'nC5', color = 'lime')
plt.plot(T, nC6/nC, label = 'nC6', color = 'pink')

plt.plot(TEvolx, nC0x/nCx, color = 'r', linestyle = ':')
plt.plot(TEvolx, nC1x/nCx, color = 'g', linestyle = ':')
plt.plot(TEvolx, nC2x/nCx, color = 'b', linestyle = ':')
plt.plot(TEvolx, nC3x/nCx, color = 'orange', linestyle = ':')
plt.plot(TEvolx, nC4x/nCx, color = 'purple', linestyle = ':')
plt.plot(TEvolx, nC5x/nCx, color = 'lime', linestyle = ':')
plt.plot(TEvolx, nC6x/nCx, color = 'pink', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)
#plt.legend()

plt.tight_layout()

plt.savefig('C_HeH_hfv.png')

plt.show()





