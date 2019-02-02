from read_molpro import read_mrci_energies as rme
from read_molpro import read_lsop_socE, LS_matrix_elements_xml, split,lsop_read_mod
import numpy as np


#au=27.2114
cminv=219474.63068
au_to_cminv= cminv

molpro_file='Example/all_r9.xml'
m = rme('Example/all_r9.xml')

#molpro_file='Example/m6_r9.xml'
#m = rme('Example/m6_r9.xml')



diag_e = np.array(m['reference Energies'])

cf_h = np.diag(diag_e-np.min(diag_e)) * au_to_cminv
#print(cf_h)

#LSOP,soc_e = read_lsop_socE(molpro_file)
LSOP,SOC_E = lsop_read_mod(molpro_file)
print(LSOP.shape)
#m6lsop1=LSOP[0: 9:1, 0: 9:1]
m6lsop=LSOP[0: 63:1, 0: 63:1]
m4lsop=LSOP[63:108:1, 63:108:1]
m2lsop=LSOP[108:135:1, 108:135:1]
m0lsop=LSOP[135:144:1, 135:144:1]
print(m6lsop.shape,m4lsop.shape,m2lsop.shape,m0lsop.shape)
#jprint(m0lsop)

#print(type(LSOP))
#print(LSOP9x9)

# LSZ=LS_matrix_elements_xml(molpro_file,'LSZ')
# LSOP_1_1_ex=LSZ+cf_h
#
# LSX=LS_matrix_elements_xml(molpro_file,'LSX')
# LSY=LS_matrix_elements_xml(molpro_file,'LSY')
#
# LSXY=LSX+LSY
# LSXZ=LSX-LSZ
#print(LSXY)


# print(LSOP.shape)
# #print(LSOP)
m6_LSOP_s=np.array(split(m6lsop,9,9))
m4_LSOP_s=np.array(split(m4lsop,9,9))
m2_LSOP_s=np.array(split(m2lsop,9,9))
m0_LSOP_s=np.array(split(m0lsop,9,9))
print(m6_LSOP_s.shape,m4_LSOP_s.shape,m2_LSOP_s.shape,m0_LSOP_s.shape)
#
# # print(np.array(LSOP_s[1]).shape)
# # print(np.all(LSOP_s==np.complex(0.0,0.0)))
# # fac=LSOP_s[0]/LSOP_s[8]
# # print(fac[0])
#
n_sub_mat=24
print(m6_LSOP_s[0].shape)
#
k=0
for j in range(len(m6_LSOP_s[0])):
    k+=1
    for i in range(len(m6_LSOP_s[0])):
        #if LSOP_s[n_sub_mat][i][j].imag != 0:
            print("%2d %2d %8.3f +( %8.3f j) %2d %2d %8.3f +( %8.3f j) %8.3f %8.3f"
                  %(i+1,j+1,
                    m6_LSOP_s[0][i][j].real,
                    m6_LSOP_s[0][i][j].imag,
                    i+1+9*n_sub_mat,
                    j+1+9*n_sub_mat,
                    m6_LSOP_s[n_sub_mat][i][j].real,
                    m6_LSOP_s[n_sub_mat][i][j].imag,
                    m6_LSOP_s[0][i][j].real/m6_LSOP_s[n_sub_mat][i][j].real,
                    m6_LSOP_s[0][i][j].imag/m6_LSOP_s[n_sub_mat][i][j].imag))
            k+=1

#LSOP_1_1  = LSOP[0][0: 9:1, 0: 9:1]
#LSOP_10_1 = LSOP[0][9:18:1, 0: 9:1]
#LSOP_10_10= LSOP[0][9:18:1, 9:18:1]

# W,V = np.linalg.eig(LSOP9x9_ex)
# print(np.sort(W-np.min(W)))

#print(LSZ)
# for i in range(9):
#     for j in range(9):
#         print('<{:2d}|LSZ|{:2d}> = {:.2f}'.format(j+1,i+1,LSZ[j,i]),
#               ' <{:2d}|LSZ|{:2d}> = {:.2f}'.format(j+1,i+1,LSOP9x9[j, i]))
#for i in range(9):
#    for j in range(9):
#        print(i+10,j+10,'{:.2f} {:.2f} {:.2f}'.format(LSOP_10_10[i, j],LSX[i, j], LSZ[i, j]))
        # if LSOP_10_1[i,j].real-LSXY[i,j].real >0.1 or LSOP_10_1[i,j].imag-LSXY[i,j].imag > 0.1:
        #
        #     print(i,j,'{:0.2f} {:0.2f}'.
        #           format(LSOP_10_1[i,j].real-LSXY[i,j].real,
        #                  LSOP_10_1[i,j].imag-LSXY[i,j].imag))