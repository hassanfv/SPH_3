
import numpy as np


#----- getAtmNum
def getAtmNum(iD):
  iDlist = np.array(['C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe'])
  AtmNumlist = [6, 7, 8, 10, 12, 14, 16, 20, 26]
  n = np.where(iD == iDlist)[0][0]
 
  return AtmNumlist[n]



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


ActiveElements = ['He', 'C', 'N', 'O', 'Ne'] # H is excluded as its abundance relative to H is 1.0 !


strx = ''

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


#-------- Constructing section for Result from "test_primordial_hdf5_v2.py" code ---------




print(strx)








