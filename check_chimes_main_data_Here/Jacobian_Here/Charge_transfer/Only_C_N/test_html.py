
import h5py
import numpy as np

# Define your lists
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

# Start HTML file
with open("reaction_data.html", "w") as html_file:
    html_file.write("<html><head><title>Chemical Reactions</title></head><body>")
    html_file.write("<h1>Chemical Reaction Data</h1>")

    # Open the HDF5 file
    with h5py.File('chimes_main_data.hdf5', 'r') as file:
        N_reactions = file['T_dependent/N_reactions'][:]
        reactants = file['T_dependent/reactants'][:]
        products = file['T_dependent/products'][:]

        for i in range(len(elmList)):
            nt = np.where(reactants[:, 0] == i)[0]   # Index for reactions involving element i
            N = len(nt)

            html_file.write(f"<h2>{elmList[i]}</h2>")

            for j in range(N):
                reac = reactants[nt[j], :]
                prod = products[nt[j], :]
                str1 = ' + '.join([elmList[k] for k in reac if k != -1]) + " ----> "
                str2 = ' + '.join([elmList[k] for k in prod if k != -1])

                label = ''
                if any(elmList[k] in Mol for k in np.concatenate((reac, prod))):
                    label = '----> MOLECULES involved !!!'
                    html_file.write(f"<p style='color:black;'>{str1}{str2} {label}</p>")
                else:
                    html_file.write(f"<p>{str1}{str2}</p>")
                    
                if any(elmList[k] in Mol for k in np.concatenate((reac, prod)) if k != -1):
                  label = '----> MOLECULES involved !!!'
                  html_file.write(f"<p style='color:red;'>{str1}{str2} {label}</p>")
                else:
                  html_file.write(f"<p style='color:black;'>{str1}{str2}</p>")

    html_file.write("</body></html>")








