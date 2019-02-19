# import numpy as np
# import numpy.matlib
# import xml.etree.ElementTree as ET
#
#
# file = 'Example/LSZ/6666.xml'
#
# # Read nstates
# f = open(file).readlines()
# for num, line in enumerate(f):
#     if 'Wavefunction restored from record' in line:
#         nstates = int(line.split()[8][-1])
#         break
# nstates=9
#
# so_mat=numpy.matlib.zeros((nstates,nstates),dtype=np.complex128)
# #so_mat=numpy.zeros((nstates,nstates),dtype=np.complex128)
# # bra = {'S': 3.0,'Ms': 2}
# # bra = {'S': 3.0,'Ms': 3}
#
# # f=fopen(file,'r')
# # for line in f :
# #     if 'Symmetry of bra wavefunction:' in line and 'S='+str(bra('S')) in line:
#
# tree = ET.parse(file)
# root=tree.getroot()
# if len(root) == 1:
#     for index,child in enumerate(root[0]):
#         if child.get('commandset')=='CIPRO':
#             ind = index
# # print(ind)
# LS_root=root[0][ind]

import numpy as np
from read_tran_ls import get_LS_mat_for_given_S as lsop_s, get_fac, LS_tran_mat_elem
from read_molpro import lsop_read_mod, split

S = 3
nstates = 9



#print(LSOP_org_split.shape)

lsop=[]

for S in [3, 2, 1, 0] :
    #for S in [3, 2, 1, 0 ] :

    if S == 0 :
        lsop_temp=np.zeros((nstates,nstates),dtype=np.complex128)
        lsop.append(lsop_temp)
    else :

        M = 2*S + 1
        org_file='Example/cas_soc_m/m'+str(2*S)+'_r'+str(nstates)+'.xml'
        LSOP_org,SOC_E_org = lsop_read_mod(org_file)

        LSOP_org_split=np.array(split(LSOP_org,nstates,nstates))

        z_prefix = str(2 * S) + str(2 * S) + str(2 * S) + str(2 * S) + '.xml'
        xy_prefix = str(2 * S) + str(2 * S) + str(2 * S - 2) + str(2 * S) + '.xml'

        file = 'Example/LSX/' + xy_prefix
        LSX = np.array(LS_tran_mat_elem(file, 'LSX', S))
        file = 'Example/LSY/' + xy_prefix
        LSY = np.array(LS_tran_mat_elem(file, 'LSY', S))
        LSXY = LSX + LSY

        file = 'Example/LSZ/' + z_prefix
        LSZ = np.array(LS_tran_mat_elem(file, 'LSZ', S))

        w3j=[]
        a=[]
        b=[]
        for i in range(0,M,1):
            a.append(get_fac(LSOP_org_split[0],LSOP_org_split[i*(M+1)]))
        for i in range(0,M-1,1):
            b.append(get_fac(LSOP_org_split[M],LSOP_org_split[M+i*(M+1)],p=False))

        w3j.append(a)
        w3j.append(b)

        lsop.append(lsop_s(S,nstates,LSX,LSY,LSZ,w3j))


#LSOP =
