
import numpy as np
import pickle


#----- getAtmNum
def getAtmNum(iD):
  iDlist = np.array(['C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe'])
  AtmNumlist = [6, 7, 8, 10, 12, 14, 16, 20, 26]
  n = np.where(iD == iDlist)[0][0]
 
  return AtmNumlist[n]


allSpeciesList = np.array([
    "elec", "H0", "H1", "Hm", "He0", "He1", "He2", "C0", "C1", "C2",
    "C3", "C4", "C5", "C6", "Cm", "N0", "N1", "N2", "N3", "N4",
    "N5", "N6", "N7", "O0", "O1", "O2", "O3", "O4", "O5", "O6",
    "O7", "O8", "Om", "Ne0", "Ne1", "Ne2", "Ne3", "Ne4", "Ne5",
    "Ne6", "Ne7", "Ne8", "Ne9", "Ne10", "Mg0", "Mg1", "Mg2", "Mg3",
    "Mg4", "Mg5", "Mg6", "Mg7", "Mg8", "Mg9", "Mg10", "Mg11", "Mg12",
    "Si0", "Si1", "Si2", "Si3", "Si4", "Si5", "Si6", "Si7", "Si8",
    "Si9", "Si10", "Si11", "Si12", "Si13", "Si14", "S0", "S1", "S2",
    "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12",
    "S13", "S14", "S15", "S16", "Ca0", "Ca1", "Ca2", "Ca3", "Ca4",
    "Ca5", "Ca6", "Ca7", "Ca8", "Ca9", "Ca10", "Ca11", "Ca12", "Ca13",
    "Ca14", "Ca15", "Ca16", "Ca17", "Ca18", "Ca19", "Ca20", "Ca21", "Fe0",
    "Fe1", "Fe2", "Fe3", "Fe4", "Fe5", "Fe6", "Fe7", "Fe8", "Fe9",
    "Fe10", "Fe11", "Fe12", "Fe13", "Fe14", "Fe15", "Fe16", "Fe17",
    "Fe18", "Fe19", "Fe20", "Fe21", "Fe22", "Fe23", "Fe24", "Fe25",
    "Fe26", "H2", "H2p", "H3p", "OH", "H2O", "C2", "O2", "HCOp", "CH",
    "CH2", "CH3p", "CO", "CHp", "CH2p", "OHp", "H2Op", "H3Op", "COp", "HOCp", "O2p"
])



# See Table_1 in Wiersma et al - 2009, 393, 99–107
elemAbund = {
    "H": 1,
    "He": 0.1,
    "C": 2.46e-4,
    "N": 8.51e-5,
    "O": 4.90e-4,
    "Ne": 1.00e-4,
    "Mg": 3.47e-5,
    "Si": 3.47e-5,
    "S": 1.86e-5,
    "Ca": 2.29e-6,
    "Fe": 2.82e-5
}


colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#1f77b4', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896',
    '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d',
    '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
    '#637939', '#8ca252'  # Add more colors if needed
]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!! Add more elements here if you want to extend your model !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ActiveElements = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe'] # H is excluded as its abundance relative to H is 1.0 !


strx = ''

strx += 'TA = time.time()\n\n'

strx += 'print("Running ...")\n\n'

strx += 'nH = 1000.0\n\n'

for elmId in ActiveElements:
  strx += f'{elmId}_solar = 10**({np.log10(elemAbund[elmId]):.2f})\n'
  strx += f'n{elmId} = {elmId}_solar * nH\n\n'


strx += 'T_i = 10**7.00\n\n'

strx += 'nHm_i = 1e-5 * nH\n'
strx += 'nH0_i = 0.001 * nH\n'
strx += 'nHp_i = nH - nH0_i\n\n'

strx += 'nHe0_i = 0.0001 * nHe\n'
strx += 'nHep_i = 0.001 * nHe\n'
strx += 'nHepp_i= nHe - nHe0_i - nHep_i\n\n'



elmList = ActiveElements[1:]

print(elmList)
print()

for elm in elmList:
    
  AtmNum = getAtmNum(elm)
  spec_list = [elm+str(i) for i in range(AtmNum+1)] # current species list ---> Contents are like C0, C1, ... C6 !!!
  
  print(spec_list)
  
  if elm == 'C':
    spec_list += ['Cm']
  if elm == 'O':
    spec_list += ['Om']

  NegIons = ['Cm', 'Om']

  for i, x in enumerate(spec_list):

    if i != AtmNum:
      initAbund = f'1e-{AtmNum - i}'
      
      if AtmNum - i < 2:
        initAbund = f'1e-2'
      
      if x in NegIons:
        initAbund = f'1e-6'
        
      strx += f'n{x}_i = {initAbund} * n{elm}\n'


  #--- adding the last transtion abundance, e.g. nC6_i = nC - (nCm_i + nC0_i + nC1_i + nC2_i + nC3_i + nC4_i + nC5_i)
  strx += f'n{spec_list[AtmNum]}_i = n{elm} - ('
  for k in range(len(spec_list)):
    if k != AtmNum:
      strx += f'n{spec_list[k]}_i + '
  strx = strx[:-3] # removing the blank and + sign
  strx += ')\n\n'


#------ Constructing the y0 list ---------
strx += 'y0 = [\n'

strx += 6*' ' + 'nH0_i, nHp_i, nHm_i, nHe0_i, nHep_i, nHepp_i, \n' + 6 * ' '

kk = 0

for elm in elmList:
    
  AtmNum = getAtmNum(elm)
  spec_list = [elm+str(i) for i in range(AtmNum+1)] # current species list ---> Contents are like C0, C1, ... C6 !!!
  
  if elm == 'C':
    spec_list += ['Cm']
  if elm == 'O':
    spec_list += ['Om']

  for x in spec_list:
    strx += f'n{x}_i, '
    kk += 1
    if not (kk % 6): # to break the lines for avoiding long lines!
      strx += '\n' + 6 * ' '

strx += '\n' + 6 * ' '
strx += 'T_i\n'
strx += 6 * ' ' + ']\n\n'
#------------ End of the y0 list ------------

strx += 'A_v = 1.0\n'
strx += 'G0 = 0.01\n'
strx += 'dust_ratio = 0.01\n\n'

strx += 't_span = (1*3.16e7, 10000*3.16e7)\n\n'

strx += 'solution = solve_ivp(func, t_span, y0, method="LSODA", dense_output=True)\n\n'

strx += 't = np.linspace(t_span[0], t_span[1], 50000) # This 10000 is not years, it is the number of points in linspace !!!!\n'
strx += 'y = solution.sol(t)\n\n'

strx += 't_yrs = t / 3.16e7\n\n'

strx += 'nH0  = y[0, :]\n'
strx += 'nHp  = y[1, :]\n'
strx += 'nHm  = y[2, :]\n'
strx += 'nHe0 = y[3, :]\n'
strx += 'nHep = y[4, :]\n'
strx += 'nHepp= y[5, :]\n\n'

jj = 6
for elm in elmList:
    
  AtmNum = getAtmNum(elm)
  spec_list = [elm+str(i) for i in range(AtmNum+1)] # current species list ---> Contents are like C0, C1, ... C6 !!!
  
  if elm == 'C':
    spec_list += ['Cm']
  if elm == 'O':
    spec_list += ['Om']

  for x in spec_list:
    strx += f'n{x} = y[{jj}, :]\n'
    jj += 1
  strx += '\n\n'

strx += f'T = y[{jj}, :]\n\n'


#-------- Constructing section for Result from CHIMES code ---------
# Note that the abundances are already multiplied by nH when the pickle file was created
strx += 'with open ("chimesData.pkl", "rb") as f:\n'
strx += '  data = pickle.load(f)\n\n'

strx += 'TEvolx = data["TempEvol"]\n'
strx += 'AbundEvol = data["chimesAbundEvol"]\n'
strx += 't_Arr_in_yrsx = data["t_Arr_in_yrs"]\n\n'

strx+= 'nH0x   = AbundEvol[1, :]\n'
strx+= 'nHpx   = AbundEvol[2, :]\n'

strx+= 'nHe0x  = AbundEvol[4, :]\n'
strx+= 'nHepx  = AbundEvol[5, :]\n'
strx+= 'nHeppx = AbundEvol[6, :]\n\n'

#---- Now other elements ----
for elm in elmList:
    
  AtmNum = getAtmNum(elm)
  spec_list = [elm+str(i) for i in range(AtmNum+1)] # current species list ---> Contents are like C0, C1, ... C6 !!!
  
  if elm == 'C':
    spec_list += ['Cm']
  if elm == 'O':
    spec_list += ['Om']

  for x in spec_list:
    ndx = np.where(allSpeciesList == x)[0][0]
    strx += f'n{x}x = AbundEvol[{ndx}, :]\n'
    
  strx += f'n{elm}x = '
  for x in spec_list: # including this: nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x
    strx += f'n{x}x + '
  strx = strx [:-2] + '\n' 
  strx += '\n'

#----- Plotting section --------
strx += 'plt.figure(figsize=(10, 5))\n'

strx += 'plt.scatter(t_yrs, np.log10(T), s=2, color="k", label="my code")\n'
strx += 'plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s=2, color="orange", label="chimes result", linestyle="--")\n'
strx += 'plt.xlim(0, 10000)\n'
strx += 'plt.ylim(1, 8)\n'
strx += 'plt.legend()\n'
strx += 'plt.savefig("Temp_vs_time.png", dpi = 300)\n'
strx += 'plt.close()\n\n'

strx += "plt.plot(t_yrs, nHe0, color = 'r', label = 'nHe0')\n"
strx += "plt.plot(t_yrs, nHep, color = 'g', label = 'nHep')\n"
strx += "plt.plot(t_yrs, nHepp, color = 'b', label = 'nHepp')\n"
strx += "plt.plot(t_Arr_in_yrsx, nHe0x, color = 'r', label = 'nHe0 - chimes', linestyle = ':')\n"
strx += "plt.plot(t_Arr_in_yrsx, nHepx, color = 'g', label = 'nHep - chimes', linestyle = ':')\n"
strx += "plt.plot(t_Arr_in_yrsx, nHeppx, color = 'b', label = 'nHepp - chimes', linestyle = ':')\n"
strx += "plt.xlim(0, 10000)\n"
strx += "plt.ylim(1e-8, 300)\n"
strx += "plt.yscale('log')\n"
strx += "plt.title('solve_ivp')\n"
strx += "plt.legend()\n"
strx += 'plt.savefig("nH_He_vs_time.png", dpi = 300)\n'
strx += 'plt.close()\n\n'


#---- Now other elements ----
for elm in elmList:
    
  AtmNum = getAtmNum(elm)
  spec_list = [elm+str(i) for i in range(AtmNum+1)] # current species list ---> Contents are like C0, C1, ... C6 !!!
  
  if elm == 'C':
    spec_list += ['Cm']
  if elm == 'O':
    spec_list += ['Om']

  for i, x in enumerate(spec_list):
    strx += f"plt.plot(T, n{x}/n{elm}, label = 'n{x}', color = '{colors[i]}')\n"
    strx += f"plt.plot(TEvolx, n{x}x/n{elm}x, label = 'n{x}x', color = '{colors[i]}', linestyle = ':')\n"
  strx += "plt.yscale('log')\n"
  strx += "plt.xscale('log')\n"
  strx += "plt.ylim(2e-3, 1.2)\n"
  strx += "plt.xlim(1e4, 1e7)\n"
  strx += f'plt.savefig("n{elm}_vs_time.png", dpi = 300)\n'
  strx += 'plt.close()\n\n'

strx += "print('Done !!!')\n\n"

strx += "print('Elapsed time = ', time.time() - TA)\n\n"

print(strx)

with open('plotHelper.py', 'w') as file:
    file.write(strx)

print('\n\n')
print('------------------------------------------------')
print()
print('        File saved to plotHelper.py file       ')
print()
print('------------------------------------------------')
print()






