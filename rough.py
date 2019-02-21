import numpy as np
import numpy.matlib
import xml.etree.ElementTree as ET


file = 'Example/LSY/ls_tran.xml'

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

tree = ET.parse(file)
root=tree.getroot()
if len(root) == 1:
    for index,child in enumerate(root[0]):
        if child.get('commandset')=='CIPRO':
            ind = index
# print(ind)
LS_root=root[0][ind]