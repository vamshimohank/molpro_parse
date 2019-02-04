import numpy as np
from read_molpro import split
from read_molpro import read_mrci_energies as rme


from read_tran_ls import *

S=1
M=2*S+1
file='Example/LSX/2202.xml'
LSX=np.array(LS_tran_mat_elem(file,'LSX',S))
file='Example/LSY/2202.xml'
LSY=np.array(LS_tran_mat_elem(file,'LSY',S))
file='Example/LSZ/2222.xml'
LSZ=np.array(LS_tran_mat_elem(file,'LSZ',S))


LSXY=LSX+LSY

from read_molpro import LS_matrix_elements_xml as LSMX
from read_molpro import lsop_read_mod

LSZ_org=LSMX('Example/m2_r9.xml','LSZ',1)
LSX_org=LSMX('Example/m2_r9.xml','LSX',1)
LSOP_org,SOC_E_org = lsop_read_mod('Example/mod_m2_r9.xml')

LSOP_org_11=LSOP_org[0:9:1,0:9:1]
LSOP_org_1010=LSOP_org[9:18:1,9:18:1]

LSOP_org_split=np.array(split(LSOP_org,9,9))

print(LSOP_org_split.shape)

mat=LSOP_org_split[8]
#mat=LSX_org
mat2=LSZ

a=[]
for i in range(0,M,1):
    print(i*M+1)
    a.append(get_fac(LSOP_org_split[0],LSOP_org_split[i*(M+1)]))

print(a)

#get_fac(LSOP_org_split[7],LSOP_org_split[5])


LSOP_ex=np.zeros((27,27),dtype=np.complex128)

m = rme('Example/m2_r9.xml')
au_to_cminv = 219474.63068
diag_e = np.array(m['reference Energies'])
cf_h = np.diag(diag_e-np.min(diag_e)) * au_to_cminv

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