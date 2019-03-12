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

file='Example/Na3Co2SbO6/soc.out'
from read_molpro import read_lsop_socE as rlsop

LSOP,E_SOC=rlsop(file)
print('LSOP.shape',LSOP.shape)
