import numpy as np
from read_molpro import split
from read_molpro import read_mrci_energies as rme

from read_molpro import LS_matrix_elements_xml as LSMX
from read_molpro import lsop_read_mod

from read_tran_ls import *

nstates = 9
S = 2
M = 2*S+1
z_prefix=str(2*S)+str(2*S)+str(2*S)+str(2*S)+'.xml'
xy_prefix=str(2*S)+str(2*S)+str(2*S-2)+str(2*S)+'.xml'

file = 'Example/LSX/'+xy_prefix
LSX = np.array(LS_tran_mat_elem(file,'LSX',S))
#file = 'Example/LSY/6646.xml'
file = 'Example/LSY/'+xy_prefix
LSY = np.array(LS_tran_mat_elem(file,'LSY',S))
file = 'Example/LSZ/'+z_prefix
#file = 'Example/LSZ/6666.xml'
LSZ = np.array(LS_tran_mat_elem(file,'LSZ',S))

print(LSZ[1][0])


LSXY = LSX+LSY

#LSZ_org = LSMX('Example/m6_r9.xml','LSZ',1)
#LSX_org = LSMX('Example/m6_r9.xml','LSX',1)
org_file='Example/m'+str(2*S)+'_r9.xml'
LSOP_org,SOC_E_org = lsop_read_mod(org_file)

LSOP_org_split=np.array(split(LSOP_org,nstates,nstates))

print(LSOP_org_split.shape)

a=[]
for i in range(0,M,1):
    a.append(get_fac(LSOP_org_split[0],LSOP_org_split[i*(M+1)]))

print(a)

k=7
b=[]
for i in range(0,M-1,1):
    b.append(get_fac(LSOP_org_split[M],LSOP_org_split[M+i*(M+1)],p=False))
print(b)

c=[]
for i in range(0,M-1,1):
    c.append(get_fac(LSOP_org_split[1],LSOP_org_split[1+i*(M+1)],p=False))
print(c)


m = rme('Example/m6_r9.xml')
au_to_cminv = 219474.63068
#diag_e = np.array(m['reference Energies'])
#cf_h = np.diag(diag_e-np.min(diag_e)) * au_to_cminv
diag_e = np.diagonal(LSOP_org[0:9, 0:9])
cf_h = np.diag(diag_e)



LSOP_ex=np.zeros((M*nstates,M*nstates),dtype=np.complex128)

LSOP_ex_split=np.array(split(LSOP_ex,nstates,nstates))

LSOP_ex_split[0] = LSZ
LSOP_ex_split[-1] = -LSZ
for l in range(M):
    if l != S :
        LSOP_ex_split[l*(M+1)] = LSZ/a[l]

for l1 in range(M-1):
    LSOP_ex_split[M+l1*(M+1)] = LSXY/b[l1]

for l2 in range(M-1):
    LSOP_ex_split[1+l2*(M+1)] = -np.conj(LSXY)/b[l1]

# from scipy.linalg import block_diag
# X=[LSOP_ex_split[i*(M+1)] for i in range(M)]
# print(len(X))
# LSOP_ex_n = block_diag(*X)
# print(LSOP_ex_n.shape)

LSOP_ex_n = np.zeros((M*nstates,M*nstates),dtype=np.complex128)
for l in range(M):
    if l == S :
        LSOP_ex_n[l*9:(l+1)*9, l*9:(l+1)*9] += cf_h
    else :
        LSOP_ex_n[l*9:(l+1)*9, l*9:(l+1)*9] = cf_h + LSZ/a[l]
#LSOP_ex_n = block_diag(LSOP_ex_split[0],LSOP_ex_split[1])

#LSOP_ex_n[ 9:18, 0: 9] = LSOP_ex_split[7]
#LSOP_ex_n[18:27, 9:18] = LSOP_ex_split[15]
#LSOP_ex_n[27:36,18:27] = LSOP_ex_split[15]

for i in range(M-1):
    #LSOP_ex_n[(i+1)*9:(i+2)*9,i*9:(i+1)*9] = LSOP_ex_split[M+i*(M+1)]
    LSOP_ex_n[(i+1)*9:(i+2)*9,i*9:(i+1)*9] = LSXY/b[i]
    LSOP_ex_n[i*9:(i+1)*9,(i+1)*9:(i+2)*9] = -np.conj(LSXY)/b[i]

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

get_diff(LSOP_ex_n,LSOP_org)

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