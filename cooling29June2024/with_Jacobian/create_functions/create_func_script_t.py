
# In this version, i.e. ***_const3.py, H and He ODEs are also created automatically!

import h5py
import numpy as np


elmList = [
            "e", "HI", "HII", "Hm", "HeI", "HeII", "HeIII", "CI", "CII", "CIII",
            "CIV", "CV", "CVI", "CVII", "Cm", "NI", "NII", "NIII", "NIV", "NV",
            "NVI", "NVII", "NVIII", "OI", "OII", "OIII", "OIV", "OV", "OVI", "OVII",
            "OVIII", "OIX", "Om", "NeI", "NeII", "NeIII", "NeIV", "NeV", "NeVI",
            "NeVII", "NeVIII", "NeIX", "NeX", "NeXI", "MgI", "MgII", "MgIII", "MgIV",
            "MgV", "MgVI", "MgVII", "MgVIII", "MgIX", "MgX", "MgXI", "MgXII", "MgXIII",
            "SiI", "SiII", "SiIII", "SiIV", "SiV", "SiVI", "SiVII", "SiVIII", "SiIX",
            "SiX", "SiXI", "SiXII", "SiXIII", "SiXIV", "SiXV", "SI", "SII", "SIII",
            "SIV", "SV", "SVI", "SVII", "SVIII", "SIX", "SX", "SXI", "SXII", "SXIII",
            "SXIV", "SXV", "SXVI", "SXVII", "CaI", "CaII", "CaIII", "CaIV", "CaV",
            "CaVI", "CaVII", "CaVIII", "CaIX", "CaX", "CaXI", "CaXII", "CaXIII", "CaXIV",
            "CaXV", "CaXVI", "CaXVII", "CaXVIII", "CaXIX", "CaXX", "CaXXI", "FeI",
            "FeII", "FeIII", "FeIV", "FeV", "FeVI", "FeVII", "FeVIII", "FeIX", "FeX",
            "FeXI", "FeXII", "FeXIII", "FeXIV", "FeXV", "FeXVI", "FeXVII", "FeXVIII",
            "FeXIX", "FeXX", "FeXXI", "FeXXII", "FeXXIII", "FeXXIV", "FeXXV", "FeXXVI",
            "FeXXVII", "H2", "H2p", "H3p", "OH", "H2O", "C2", "O2", "HCOp", "CH",
            "CH2", "CH3p", "CO", "CHp", "CH2p", "OHp", "H2Op", "H3Op", "COp", "HOCp", "O2p"
          ]


Mol = ["H2", "H2p", "H3p", "OH", "H2O", "C2", "O2", "HCOp", "CH",
       "CH2", "CH3p", "CO", "CHp", "CH2p", "OHp", "H2Op", "H3Op", "COp", "HOCp", "O2p"]


roman_num = [
                  "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                  "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",
                  "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII"
                 ]


grain_rec_reac = np.array(['HII', 'HeII', 'CII', 'OII', 'SiII', 'FeII', 'MgII', 'SII', 'CaII'])#, 'CaIII']) # for CIII_to_CaII I hard-coded it see below!!
grain_rec_prod = np.array(['HI', 'HeI', 'CI', 'OI', 'SiI', 'FeI', 'MgI', 'SI', 'CaI'])#, 'CaII'])

#----- constant rates ------- !!!!!!!!!!! DO NOT FORGET TO HARD-CODE Cm and Om !!!!!!!!!!!!!!!!
constList = [["CII", "SiI", "CI", "SiII"],
             ["OI", "e", "Om", ""],
             ["CII", "MgI", "CI", "MgII"],
             ["NII", "MgI", "NI", "MgII"],
             ["CI", "e", "Cm", ""],
             ["SII", "MgI", "SI", "MgII"],
             ["FeI", "SiII", "FeII", "SiI"],
             ["FeI", "CII", "FeII", "CI"],
             ["SiII", "MgI", "SiI", "MgII"]]




#----- getAtmNum ---> get Atomic Number from elements name or iD!
def getAtmNum(iD):
  iDlist = np.array(['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe'])
  AtmNumlist = [1, 2, 6, 7, 8, 10, 12, 14, 16, 20, 26]
  n = np.where(iD == iDlist)[0][0]
 
  return AtmNumlist[n]
  


res1_list = []
res2_list = []

# Open the HDF5 file
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  N_reactions = file['T_dependent/N_reactions'][:]
  print('N_reactions = ', N_reactions)
  print()
  
  reactants = file['T_dependent/reactants'][:]
  print("reactants.shape':", reactants.shape)
  products = file['T_dependent/products'][:]
  print("products.shape':", reactants.shape)
  print()

  const_rates = file['constant/rates'][:] # Seems they are just charge exchange and no cooling actually happen as a result of these processes!

  reactants_grain_recomb = file['grain_recombination/reactants'][:]
  products_grain_recomb = file['grain_recombination/products'][:]


#---------- 
elm = 'H'   
AtmNum = getAtmNum(elm)
spec_list = [elm+roman_num[i] for i in range(AtmNum+1)]
spec_iState = [i for i in range(AtmNum+1)] # current species ionization states
spec_list += ['Hm']
spec_iState += [-1]
#---------- 
elm = 'He'   
AtmNum = getAtmNum(elm)
spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
spec_list += spec_list2
spec_iState += spec_iState2
#---------- 
elm = 'C'   
AtmNum = getAtmNum(elm)
spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
spec_list += spec_list2
spec_iState += spec_iState2
spec_list += ['Cm']
spec_iState += [-1]
#---------- 

elm = 'N'   
AtmNum = getAtmNum(elm)
spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
spec_list += spec_list2
spec_iState += spec_iState2

#---------- 
elm = 'O'   
AtmNum = getAtmNum(elm)
spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
spec_list += spec_list2
spec_iState += spec_iState2
spec_list += ['Om']
spec_iState += [-1]
#----------------------------------

if False:
  #---------- 
  elm = 'Ne'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------

  #---------- 
  elm = 'Mg'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------

  #---------- 
  elm = 'Si'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------

  #---------- 
  elm = 'S'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------

  #---------- 
  elm = 'Ca'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------

  #---------- 
  elm = 'Fe'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)]
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------

print()
print('-----------> spec_list <----------')
print(spec_list)
print()
print('-----------> spec_iState <----------')
print(spec_iState)


# Writing to text files
file1 = open('func.py', 'w')

tmp = '#----- func\n'
file1.write(tmp)
tmp = 'def func(t, y):\n\n  '
file1.write(tmp)

#----------- here we create: #nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, nCV, nCVI, nCVII, nCm, T = y
tmp = ''
kk = 1
for x in spec_list:
  tmp += f'n{x}, '
  if not (kk % 10):
    tmp += "\\"
    tmp += f'\n{2 * " "}'
  kk += 1
tmp += 'T = y\n\n'
file1.write(tmp)

#----------------
tmp = '  Tx = np.log10(T)\n\n'
file1.write(tmp)

#---- nelec or ne section ----
str_ne = f'  ne = (\n{6 * " "}' # f'  ne = (\n{8 * " "} 1 * nHII - nHm + (nHeII + 2.0 * nHeIII)'
kk = 1
for k, x in zip(spec_iState, spec_list):
  if k !=0:
    str_ne += f' + {k} * n{x}'
    if not (kk % 6):
      str_ne += f'\n{6 * " "}'
    kk += 1
str_ne += f'\n{7 * " "})\n\n'
file1.write(str_ne)
#------ End of nelec section -----


#-------- ntot section -----------
tmp = f'  ntot = (\n{10 * " "} ne'
kk = 1
for x in spec_list:
  tmp += f' + n{x}'
  if not (kk % 6):
    tmp += f'\n{8 * " "}'
  kk += 1
tmp += f'\n{9 * " "})\n\n'
file1.write(tmp)
#---------------------------------

#------- grain_recombination rates ---------
for x, a in zip(reactants_grain_recomb[:, :], products_grain_recomb):
  if (elmList[x[0]] in spec_list) and (elmList[a] in spec_list):
    tmp = f'  grain_rec_{elmList[x[0]]}_to_{elmList[a]} = grain_recomb_rate("{elmList[x[0]]}", T, ne, G0, A_v, Temp, Psi)\n'
    file1.write(tmp)
#-------------------------------------------

#------- Writing the constant rates to the file ---------
j = 0
for tmpArr in constList:
  reacz = tmpArr[:2]
  prodz = tmpArr[2:]
  
  if (reacz[0] in spec_list) and ((reacz[1] in spec_list) or (reacz[1] == 'e')):
    tmp = f'\n  const_{reacz[0]}_{reacz[1]}_to_{prodz[0]}_{prodz[1]} = {const_rates[j]:.4E} # constant/rates'
    file1.write(tmp)
  j+= 1
tmp = '\n\n'
file1.write(tmp)
#--------------------------------------------------------

print()
print(f'spec_list = {spec_list}')
print()

#----------- writing the ODEs to the file ---------------
oneTimerCa = 0
for jj in range(len(spec_list)):

  iD = spec_list[jj]
  ode = f'  dn{iD}_dt = '
  Nspace = len(ode)
  checker = 0
  oneTimerConst = 0

  file1.write(ode)

  oneTimer = 0 # To make sure the grain_rec if statements are activated only once!

  for i in range(len(elmList)):

    nt = np.where(reactants[:, 0] == i)[0] 
    
    
    #print(reactants.shape)
    #print(reactants)
    print(f'elmList[i] = {elmList[i]}')
    print(f'iD = {iD}')
    print(f'nt = {nt} ---> where reactants[:, 0] == {elmList[i]}')
    print(f'reactants[nt] = {reactants[nt]}')
    print()
    print()

    N = len(nt)

    for j in range(N):
    
      a = b = c = x = y = z = ''

      reac = reactants[nt[j], :] # shape ---> 1 * 3

      #===== The reactants section ======
      a = elmList[reac[0]]
      b = elmList[reac[1]]  
      if reac[2] != -1:
        c = elmList[reac[2]]
      
      str1 = f'{a} + {b}'
      if reac[2] != -1:
        str1 = str1 + f' + {c}'

      str1 = str1 + f' ----> '

      #===== The products section ======
      prod = products[nt[j], :]

      x = elmList[prod[0]]
      if prod[1] != -1:
        y = elmList[prod[1]]
      if prod[2] != -1:
        z = elmList[prod[2]]

      str2 = f'{x}'
      if prod[1] != -1:
        str2 = str2 + f' + {y}'
      if prod[2] != -1:
        str2 = str2 + f' + {z}'
      
      if str2 == f'{x}':
        str2 = str2 + f' + γ'
      
      strx = str1 + str2
      
      #spaces = ' ' * (31 - len(strx))
      
      NNN = 31 - len(strx)
      strx = strx + NNN * ' ' + f'   (ndx = {nt[j]})'
      
      label = ''
      
      if 'γ' in strx:
        label = ' ----> radiative recombination'
      if 'e + e' in strx:
        label = ' ----> collisional ionization'
      
      if (a in Mol) or (b in Mol) or (c in Mol) or (x in Mol) or (y in Mol) or (z in Mol):
        label = '----> MOLECULES involved !!!'
      
      strz = strx + label
      
      pm = 0.0
      if label != '----> MOLECULES involved !!!':
      
        if (iD in [a, b, c]) and (a in spec_list) and ((b in spec_list) or (b == 'e')): # list of reactants.
          if checker == 0: # checker is used for the paranthesis in the begining of the line!!!! Good idea!
            tmp = f'(\n{Nspace * " "} - 10**R_{a}_to_{x}_via_{b}(Tx) * n{a} * n{b}'
            ode = ode + tmp
            checker = 1
            file1.write(tmp)
          else:
            tmp = f'\n{Nspace * " "} - 10**R_{a}_to_{x}_via_{b}(Tx) * n{a} * n{b}'
            ode = ode + tmp
            checker = 1
            file1.write(tmp)
        if (iD in grain_rec_reac) and (checker == 1) and (oneTimer == 0): # if it is the reactant then we lose it so it needs to be subtracted!
            nn = np.where(grain_rec_reac == iD)[0][0]
            tmp = f'\n{Nspace * " "} - 10**grain_rec_{iD}_to_{grain_rec_prod[nn]} * n{iD} * ne # grain_recombination'
            ode = ode + tmp
            oneTimer = 1
            file1.write(tmp)
          
        if (iD in grain_rec_prod) and (checker == 1) and (oneTimer == 0): #if it is the product then we gain it so it needs to be dded!
            nn = np.where(grain_rec_prod == iD)[0][0]
            tmp = f'\n{Nspace * " "} + 10**grain_rec_{grain_rec_reac[nn]}_to_{iD} * n{grain_rec_reac[nn]} * ne # grain_recombination'
            ode = ode + tmp
            oneTimer = 1
            file1.write(tmp)
        
        
        #--- Due to the complexiy I prefer to hard-code these two lines!!!
        if (iD == 'CaII') and (checker == 1) and (oneTimerCa == 0): # checker is needed so that it is not added as the first line because it hasn't paranthesis (!!
          tmp = f'\n{Nspace * " "} + 10**grain_rec_CaIII_to_CaII * nCaIII * ne # grain_recombination'
          ode = ode + tmp
          oneTimerCa = 1
          file1.write(tmp)
          
        if (iD == 'CaIII') and (checker == 1) and (oneTimerCa == 1): # checker is needed so that it is not added as the first line because it hasn't paranthesis (!!
          tmp = f'\n{Nspace * " "} - 10**grain_rec_CaIII_to_CaII * nCaIII * ne # grain_recombination'
          ode = ode + tmp
          oneTimerCa = 2
          file1.write(tmp)
        #-----------------------------------------------------------------
        
        if (iD in [x, y, z]) and (x in spec_list) and ((y in spec_list) or (y == 'e') or (len(y) == 0)): # list of products.
          if checker == 0:
            tmp = f'(\n{Nspace * " "} + 10**R_{a}_to_{x}_via_{b}(Tx) * n{a} * n{b}'
            ode = ode + tmp
            checker = 1
            file1.write(tmp)
          else:
            tmp = f'\n{Nspace * " "} + 10**R_{a}_to_{x}_via_{b}(Tx) * n{a} * n{b}'
            ode = ode + tmp
            checker = 1
            file1.write(tmp)
  
  s()

  #------- HII & HeII Case A recombination --------
  if iD == 'HI':
    tmp = f'\n{Nspace * " "} + 10**R_HII_to_HI_via_e_caseA(Tx) * nHII * ne # H CaseA'
    ode += tmp
    file1.write(tmp)
  if iD == 'HII':
    tmp = f'\n{Nspace * " "} - 10**R_HII_to_HI_via_e_caseA(Tx) * nHII * ne # H CaseA'
    ode += tmp
    file1.write(tmp)
  if iD == 'HeI':
    tmp = f'\n{Nspace * " "} + 10**R_HeII_to_HeI_via_e_caseA(Tx) * nHeII * ne # He CaseA'
    ode += tmp
    file1.write(tmp)
  if iD == 'HeII':
    tmp = f'\n{Nspace * " "} - 10**R_HeII_to_HeI_via_e_caseA(Tx) * nHeII * ne # He CaseA'
    ode += tmp
    file1.write(tmp)
  #-------------------------------------------------
  
  #-------------- constant/rates -------------------
  for tmpArr in constList:
    reacz = tmpArr[:2]
    prodz = tmpArr[2:]
    
    if (reacz[0] in spec_list) and ((reacz[1] in spec_list) or (reacz[1] == 'e')):
      if iD == reacz[0]:
        tmp = f'\n{Nspace * " "} - const_{reacz[0]}_{reacz[1]}_to_{prodz[0]}_{prodz[1]} * n{reacz[0]} * n{reacz[1]} # constant rate'
        ode += tmp
        file1.write(tmp)
        
      if iD == reacz[1]:
        tmp = f'\n{Nspace * " "} - const_{reacz[0]}_{reacz[1]}_to_{prodz[0]}_{prodz[1]} * n{reacz[0]} * n{reacz[1]} # constant rate'
        ode += tmp
        file1.write(tmp)

      if iD == prodz[0]:
        tmp = f'\n{Nspace * " "} + const_{reacz[0]}_{reacz[1]}_to_{prodz[0]}_{prodz[1]} * n{reacz[0]} * n{reacz[1]} # constant rate'
        ode += tmp
        file1.write(tmp)

      if iD == prodz[1]:
        tmp = f'\n{Nspace * " "} + const_{reacz[0]}_{reacz[1]}_to_{prodz[0]}_{prodz[1]} * n{reacz[0]} * n{reacz[1]} # constant rate'
        ode += tmp
        file1.write(tmp)
              
        
        
  file1.write(f'\n{Nspace * " "})\n\n')

  ode = ode + f'\n{Nspace * " "})\n\n' # this is just for printing in the terminal! 

  print()
  print(ode)
  print()
  print()


#------ Adding the Lamb line ---------
tmp = f'  Lamb = Lambda(\n{16 * " "}T, '
kk = 1
for x in spec_list:
  tmp += f'n{x}, '
  if not (kk % 10):
    tmp += f'\n{16 * " "}'
  kk += 1
tmp = tmp[:-2]
tmp += ')\n\n'
file1.write(tmp)
#--------------------------------------

#---- dne_dt section ----
str_ne = f'  dne_dt = (\n{10 * " "}' # f'  ne = (\n{8 * " "} 1 * nHII - nHm + (nHeII + 2.0 * nHeIII)'
kk = 1
for k, x in zip(spec_iState, spec_list):
  if k !=0:
    str_ne += f' + {k} * dn{x}_dt'
    if not (kk % 6):
      str_ne += f'\n{10 * " "}'
    kk += 1
str_ne += f'\n{11 * " "})\n\n'
file1.write(str_ne)
#------ End of nelec section -----

#-------- dntot_dt section -----------
tmp = f'  dntot_dt = (\n{15 * " "} dne_dt'
kk = 1
for x in spec_list:
  tmp += f' + dn{x}_dt'
  if not (kk % 6):
    tmp += f'\n{13 * " "}'
  kk += 1
tmp += f'\n{13 * " "})\n\n'
file1.write(tmp)
#---------------------------------

tmp = '  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * (Lamb + 1. / (gamma - 1.) * kB * T * dntot_dt)\n\n'
file1.write(tmp)

#-------- return section -----------
tmp = f'  return [\n{8 * " "} '
kk = 1
for x in spec_list:
  tmp += f' dn{x}_dt,'
  if not (kk % 6):
    tmp += f'\n{9 * " "}'
  kk += 1
tmp += ' dT_dt'
tmp += f'\n{9 * " "}]\n\n'
file1.write(tmp)
#---------------------------------

file1.close() # The output file name is -----> ode1.py

print()
print()
print('**********************************************************************')
print('                     !!! IMPORTANT IMPORTANT !!!')
print('**********************************************************************')
print()
print('REMOVE THIS ---> - 10**R_Hm_to_HI_via_HI(Tx) * nHm * nHI <--- FROM dnHI_dt')
print('Reason: This raction, i.e. Hm + HI ---> HI + HI + e adds HI but ')
print('my script also subtracts it because HI is present on both sides!')
print()
print('**********************************************************************')
print()
print()
print('               Result written in ---> func.py')
print()


