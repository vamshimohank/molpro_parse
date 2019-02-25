import numpy as np
from read_molpro import split
from read_molpro import read_mrci_energies as rme

from read_molpro import read_lsop_socE as r_lsop_E
from read_molpro import lsop_read_mod

from read_tran_ls import *

nstates = 9
S1 = 3
S2 = 2
M1 = 2*S1+1
M2 = 2*S2+1

# z_prefix=str(2*S)+str(2*S)+str(2*S)+str(2*S)+'.xml'
z_prefix=str(2*S1)+str(2*S2)+str(2*S2)+str(2*S2)+'_LSZ.xml'
x_prefix=str(2*S1)+str(2*S2)+str(2*S1)+str(2*S2)+'_LSX.xml'
y_prefix=str(2*S1)+str(2*S2)+str(2*S1)+str(2*S2)+'_LSY.xml'
#
file = 'Example/LSX/'+x_prefix
LSX = np.array(LS_tran_mat_elem(file,'LSX',S1,S2))
print('LSX.shape',LSX.shape)
file = 'Example/LSY/'+y_prefix
LSY = np.array(LS_tran_mat_elem(file,'LSY',S1,S2))
print('LSY.shape',LSY.shape)
file = 'Example/LSZ/'+z_prefix
LSZ = np.array(LS_tran_mat_elem(file,'LSZ',S1,S2))
print('LSZ.shape',LSZ.shape)
# #file = 'Example/LSY/6646.xml'
# file = 'Example/LSY/'+xy_prefix
# LSY = np.array(LS_tran_mat_elem(file,'LSY',S))
# file = 'Example/LSZ/'+z_prefix
# #file = 'Example/LSZ/6666.xml'
# LSZ = np.array(LS_tran_mat_elem(file,'LSZ',S))
#
# print(LSZ[1][0])


LSXY = LSX+LSY
print('LSXY.shape',LSXY.shape)

#LSZ_org = LSMX('Example/m6_r9.xml','LSZ',1)
#LSX_org = LSMX('Example/m6_r9.xml','LSX',1)

# f='Example/Na3Co2SbO6/soc.out'
f='Example/Na2Co2TeO6/cas10_soc.xml'

org_file='Example/all_r9.xml'
LSZ_org = np.array(LS_tran_mat_elem(org_file,'LSZ',S1,S2))
LSX_org = np.array(LS_tran_mat_elem(org_file,'LSX',S1,S2))
LSY_org = np.array(LS_tran_mat_elem(org_file,'LSY',S1,S2))
print("LSZ_org.shape =",LSZ_org.shape)
print("LSX_org.shape =",LSX_org.shape)
print("LSY_org.shape =",LSY_org.shape)

LSXY_org=LSX_org+LSY_org

# print_mat(LSZ_org)

# LSOP_org,SOC_E_org = lsop_read_mod(org_file)
LSOP_org,SOC_E_org = r_lsop_E(f)
LSOP_org=LSOP_org[0]
print('LSOP_org.shape',LSOP_org.shape)
print(M1*nstates,nstates*M2)

ls33=LSOP_org[0:63:1, 0:63:1]

# print_mat(ls33)
ls32=LSOP_org[0:63:1, 90:135:1]
# print_mat(ls32)
ls32_ct = LSOP_org[90:135:1, 0:63:1]
get_diff(ls32,np.transpose(ls32_ct))

#ls32_split=np.array(split(ls32,nstates,nstates)) # it looks like this is not working
#print(ls32_split.shape)


temp=-(LSX+LSY)
temp1=-np.conj(LSX+LSY)

#get_diff(ls32_ct[0:9:1, 0:9:1],LSX+LSY)
# get_diff(ls32_ct[0:9:1, 0:9:1],temp1)
#get_diff(ls32[0:9:1, 0:9:1],temp)
# temp1 = -np.conj(temp)
#temp1 = -np.conj(temp1)
# get_diff_3(LSX,np.transpose(LSY),temp1)
# get_diff_3(LSY,LSX,temp1)
# get_diff(LSXY,ls32_split[8])
# print_mat(ls32,ls32_ct)

# get_fac(ls32,ls32_ct)
# print(temp.shape,temp1.shape)
# for i in range(temp.shape[0]):
#     for j in range(temp.shape[1]):
#         print(i+63,j,temp[i][j],temp1[j][i])
# LSOP_org_split=np.array(split(LSOP_org,nstates,nstates))
#
# print(LSOP_org_split.shape)
#
# a=[]
# for i in range(0,M,1):
#     a.append(get_fac(LSOP_org_split[0],LSOP_org_split[i*(M+1)]))
#
# print(a)
#
# k=7
# b=[]
# for i in range(0,M-1,1):
#     b.append(get_fac(LSOP_org_split[M],LSOP_org_split[M+i*(M+1)],p=False))
# print(b)
#
# c=[]
# for i in range(0,M-1,1):
#     c.append(get_fac(LSOP_org_split[1],LSOP_org_split[1+i*(M+1)],p=False))
# print(c)
#
#
# m = rme('Example/m6_r9.xml')
# au_to_cminv = 219474.63068
# #diag_e = np.array(m['reference Energies'])
# #cf_h = np.diag(diag_e-np.min(diag_e)) * au_to_cminv
# diag_e = np.diagonal(LSOP_org[0:9, 0:9])
# cf_h = np.diag(diag_e)
#
#
#
# LSOP_ex=np.zeros((M*nstates,M*nstates),dtype=np.complex128)
#
# LSOP_ex_split=np.array(split(LSOP_ex,nstates,nstates))
#
# LSOP_ex_split[0] = LSZ
# LSOP_ex_split[-1] = -LSZ
# for l in range(M):
#     if l != S :
#         LSOP_ex_split[l*(M+1)] = LSZ/a[l]
#
# for l1 in range(M-1):
#     LSOP_ex_split[M+l1*(M+1)] = LSXY/b[l1]
#
# for l2 in range(M-1):
#     LSOP_ex_split[1+l2*(M+1)] = -np.conj(LSXY)/b[l1]
#
# # from scipy.linalg import block_diag
# # X=[LSOP_ex_split[i*(M+1)] for i in range(M)]
# # print(len(X))
# # LSOP_ex_n = block_diag(*X)
# # print(LSOP_ex_n.shape)
#
# LSOP_ex_n = np.zeros((M*nstates,M*nstates),dtype=np.complex128)
# for l in range(M):
#     if l == S :
#         LSOP_ex_n[l*9:(l+1)*9, l*9:(l+1)*9] += cf_h
#     else :
#         LSOP_ex_n[l*9:(l+1)*9, l*9:(l+1)*9] = cf_h + LSZ/a[l]
# #LSOP_ex_n = block_diag(LSOP_ex_split[0],LSOP_ex_split[1])
#
# #LSOP_ex_n[ 9:18, 0: 9] = LSOP_ex_split[7]
# #LSOP_ex_n[18:27, 9:18] = LSOP_ex_split[15]
# #LSOP_ex_n[27:36,18:27] = LSOP_ex_split[15]
#
# for i in range(M-1):
#     #LSOP_ex_n[(i+1)*9:(i+2)*9,i*9:(i+1)*9] = LSOP_ex_split[M+i*(M+1)]
#     LSOP_ex_n[(i+1)*9:(i+2)*9,i*9:(i+1)*9] = LSXY/b[i]
#     LSOP_ex_n[i*9:(i+1)*9,(i+1)*9:(i+2)*9] = -np.conj(LSXY)/b[i]

#print(LSOP_ex_split[8]-LSOP_ex_n[9:18, 9:18])


# ev_to_cminv = 8065.5
# e_org=np.linalg.eigvals(LSOP_org)
# e_org = e_org-np.min(e_org)
# print(e_org)
# e_ex=np.linalg.eigvals(LSOP_ex_n)
# e_ex=e_ex-np.min(e_ex)
# print(e_ex)
# Ediff=e_ex.real-e_org.real
# print(Ediff/ev_to_cminv)

# get_diff(LSOP_ex_n,LSOP_org)

"""
LSOP_ex[ 0: 9:1, 0: 9:1] =  LSZ             + cf_h      #[0]
LSOP_ex[ 9:18:1, 9:18:1] +=                 + cf_h      #[4]
LSOP_ex[18:27:1,18:27:1] = -LSZ             + cf_h      #[8]
LSOP_ex[ 9:18:1, 0: 9:1] =  LSX+LSY                     #[3]
LSOP_ex[ 0: 9:1, 9:18:1] = -np.conj(LSX+LSY)            #[1]
LSOP_ex[18:27:1, 9:18:1] =  LSX+LSY                     #[7]
LSOP_ex[ 9:18:1,18:27:1] = -np.conj(LSX+LSY)            #[5]

LSOP_ex_split=np.array(split(LSOP_ex,9,9))


# for k in range(9):
# l = k + 9
# LSOP_ex[k:l:1, k:l:1] = LSZ + cf_h
#
ev_to_cminv = 8065.5
#print(diag_e)

#print(cf_h)
#print_mat(LSOP_ex)
#get_fac(LSOP_ex,LSOP_org)
#get_diff(LSOP_ex,LSOP_org)
#print(np.diag(LSOP_org)-np.diag(LSOP_ex))

e_org = np.linalg.eigvals(LSOP_org)
e_org = e_org-np.min(e_org)
#print(e_org)

e_ex=np.linalg.eigvals(LSOP_ex)
e_ex=e_ex-np.min(e_ex)
Ediff=e_ex.real-e_org.real
print(Ediff/ev_to_cminv)

"""