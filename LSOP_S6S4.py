import numpy as np
from read_molpro import split
from read_molpro import read_mrci_energies as rme

from read_molpro import read_lsop_socE as r_lsop_E
from read_molpro import lsop_read_mod

from read_tran_ls import *

from read_tran_ls import LS_tran_mat_elem as read_tran_mat_ele

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
file = 'Example/Na2Co2TeO6/LSX/'+x_prefix
# print(file)
LSX = np.array(LS_tran_mat_elem(file,'LSX',S1,S2))
print('LSX.shape',LSX.shape)
# print_mat(LSX)
file = 'Example/Na2Co2TeO6/LSY/'+y_prefix
LSY = np.array(LS_tran_mat_elem(file,'LSY',S1,S2))
print('LSY.shape',LSY.shape)
file = 'Example/Na2Co2TeO6/LSZ/'+z_prefix
print(file)
LSZ = np.array(read_tran_mat_ele(file,'LSZ',S1,S2))
print('LSZ.shape',LSZ.shape)


LSXY = LSX+LSY
print('LSXY.shape',LSXY.shape)

f='Example/Na2Co2TeO6/cas10_soc.xml'

LSOP_org,SOC_E_org = r_lsop_E(f)
LSOP_org=LSOP_org[0]
# print('LSOP_org.shape',LSOP_org.shape)
# print(M1*nstates,nstates*M2)

ls33=LSOP_org[0:63:1, 0:63:1]

# print_mat(ls33)


ls32_ct = LSOP_org[63:108:1, 0:63:1]
ls32=LSOP_org[0:63:1, 63:108:1]
get_diff(ls32,np.transpose(np.conjugate(ls32_ct)),op='-')
print('ls32.shape=',ls32.shape)

ls32_split=np.array(blockshaped(ls32,9,9))


# temp_fac= get_fac(ls32_split[0],ls32_split[6],p=False)
# print(temp_fac)

new_ls32=np.zeros((M1*nstates,M2*nstates),dtype=np.complex128)
print(new_ls32.shape)
new_ls32_split=np.array(blockshaped(new_ls32,nstates,nstates))
#
b=[]
for i in range(0,M2,1):
    fac=get_fac(ls32_split[0],ls32_split[i*(M2+1)])
    b.append(fac)
    new_ls32_split[i*(M2+1)]=LSXY/b[i]
# print(b)

a=[]
for i in range(0,M2,1):
    fac=get_fac(ls32_split[M2],ls32_split[M2+i*(M2+1)])
    a.append(fac)
    new_ls32_split[M2+i*(M2+1)]=LSZ/a[i]
# print(a)

bs=[]
for i in range(0,M2,1):
    fac=get_fac(ls32_split[0],np.conjugate(ls32_split[2*M2+i*(M2+1)]))
    bs.append(fac)
    new_ls32_split[2*M2+i*(M2+1)]=np.conjugate(LSXY)/bs[i]
# print(bs)
# get_fac(ls32_split[0],ls32_split[34],p=True)


for i in range(M2):
    # print(new_ls32[i*9:i*9+9, i*9:i*9+9].shape,new_ls32_split[i*(M2+1)].shape)
    new_ls32[i*9:i*9+9, i*9:i*9+9] = new_ls32_split[i*(M2+1)]
    new_ls32[(i+1)*9:(i+1)*9+9, i*9:i*9+9] = new_ls32_split[M2+i*(M2+1)]
    new_ls32[(i+2)*9:(i+2)*9+9, i*9:i*9+9] = new_ls32_split[2*M2+i*(M2+1)]

new_ls32_ct=np.transpose(np.conjugate(new_ls32))
# get_diff(new_ls32,np.transpose(np.conjugate(new_ls32_ct)),op='-',p=True)

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