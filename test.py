from read_molpro import read_mrci_energies as rme
from read_molpro import read_lsop_socE, LS_matrix_elements_xml
import numpy as np


#au=27.2114
cminv=219474.63068
au_to_cminv= cminv

molpro_file='Example/m7_soc.xml'

m = rme('Example/m7_SOC.xml')

diag_e = np.array(m['reference Energies'])

cf_h = np.diag(diag_e-np.min(diag_e)) * au_to_cminv
#print(cf_h)

LSOP,soc_e = read_lsop_socE(molpro_file)

#print(type(LSOP))
#print(LSOP9x9)

LSZ=LS_matrix_elements_xml(molpro_file,'LSZ')
LSOP_1_1_ex=LSZ+cf_h

LSX=LS_matrix_elements_xml(molpro_file,'LSX')
LSY=LS_matrix_elements_xml(molpro_file,'LSY')

LSXY=LSX+LSY
LSXZ=LSX-LSZ
#print(LSXY)

LSOP_1_1  = LSOP[0][0: 9:1, 0: 9:1]
LSOP_10_1 = LSOP[0][9:18:1, 0: 9:1]
LSOP_10_10= LSOP[0][9:18:1, 9:18:1]

# W,V = np.linalg.eig(LSOP9x9_ex)
# print(np.sort(W-np.min(W)))

#print(LSZ)
# for i in range(9):
#     for j in range(9):
#         print('<{:2d}|LSZ|{:2d}> = {:.2f}'.format(j+1,i+1,LSZ[j,i]),
#               ' <{:2d}|LSZ|{:2d}> = {:.2f}'.format(j+1,i+1,LSOP9x9[j, i]))
for i in range(9):
    for j in range(9):
        print(i+10,j+10,'{:.2f} {:.2f} {:.2f}'.format(LSOP_10_10[i, j],LSX[i, j], LSZ[i, j]))
        # if LSOP_10_1[i,j].real-LSXY[i,j].real >0.1 or LSOP_10_1[i,j].imag-LSXY[i,j].imag > 0.1:
        #
        #     print(i,j,'{:0.2f} {:0.2f}'.
        #           format(LSOP_10_1[i,j].real-LSXY[i,j].real,
        #                  LSOP_10_1[i,j].imag-LSXY[i,j].imag))