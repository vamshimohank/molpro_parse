import numpy as np
import numpy.matlib
import xml.etree.ElementTree as ET


file = 'Example/all_r9.xml'

# Read nstates
# f = open(file).readlines()
# for num, line in enumerate(f):
#     if 'Wavefunction restored from record' in line:
#         nstates = int(line.split()[8][-1])
#         break
nstates=9

so_mat=numpy.matlib.zeros((nstates,nstates),dtype=np.complex128)
#so_mat=numpy.zeros((nstates,nstates),dtype=np.complex128)
# bra = {'S': 3.0,'Ms': 2}
# bra = {'S': 3.0,'Ms': 3}

# f=fopen(file,'r')
# for line in f :
#     if 'Symmetry of bra wavefunction:' in line and 'S='+str(bra('S')) in line:

# tree = ET.parse(file)
# root=tree.getroot()
# if len(root) == 1:
#     for index,child in enumerate(root[0]):
#         if child.get('commandset')=='CIPRO':
#             ind = index
# print(ind)
# LS_root=root[0][ind]
# for i in range(len(LS_root)):
#     if LS_root[i].get('name') == 'LSX' :
#         print(LS_root[i].get('name') == 'LSX')

file='Example/Na2Co2TeO6/cas10_soc.xml'
from read_molpro import read_lsop_socE as r_lsop_E

from read_tran_ls import *

LSOP,E_SOC=r_lsop_E(file)
print('LSOP.shape',LSOP.shape)

LSOP_org, SOC_E_org = r_lsop_E(file)
LSOP_org = LSOP_org[0]

# ls32_ct = LSOP_org[63:108:1, 0:63:1]
ls32 = LSOP_org[0:63:1, 63:108:1]
# get_diff(ls32, np.transpose(np.conjugate(ls32_ct)), op='-')
# print('ls32.shape=', ls32.shape)

ls32_split = np.array(blockshaped(ls32, 9, 9))
LSOP_new=build_s_sp1_lsmat(LS_org=ls32_split,S1=3,S2=2,nstates=9,files_path='Example/Na2Co2TeO6',x_prefix='6464_LSX.xml',
                           y_prefix='6464_LSY.xml',z_prefix='6444_LSZ.xml')

ls21 = LSOP_org[63:108:1, 108:135:1]
print(ls21.shape)
LSOP_new=build_s_sp1_lsmat(S1=1,S2=2,nstates=9,files_path='Example/Na2Co2TeO6',x_prefix='2424_LSX.xml',
                           y_prefix='2424_LSY.xml',z_prefix='2422_LSZ.xml')


"""
import numpy as np

array = np.array([
    [1, 1, 2, 2],
    [3, 3, 4, 4],
    [5, 5, 6, 6],
    [7, 7, 8, 8],
    [8, 8, 10, 10],
    [9, 9, 11, 11]])
print(array)

from read_molpro import split

x = np.array(blockshaped(array, 2, 2))
print(x.shape)
print(x)
"""

