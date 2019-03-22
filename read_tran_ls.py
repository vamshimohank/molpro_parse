

import numpy as np
def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def LS_tran_mat_elem(file,lstype,S1,S2):


    au=27.2114
    cminv=8.07
    autocminv= au * cminv * 1000

    au_to_cminv = 219474.6308408

    from numpy import matrixlib
    import numpy as np
    import numpy.matlib
    import xml.etree.ElementTree as ET


    # Read nstates
    f = open(file).readlines()
    for num, line in enumerate(f):
        if 'Wavefunction restored from record' in line:
            nstates = int(line.split()[8][-1])
            break
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

    method = 'MRCI'
    index = 0
    # method = 'MRCI LS_I-I(FV)'

    #print(LS_root[0].attrib)
    # print(len(LS_root))

    if 'ECP' in lstype :
        print('In ECP')
        for i in range(len(LS_root)):
            if LS_root[i].get('name') == lstype and LS_root[i].get('method') == 'MRCI LS_I-I(ECP)':

                print('<%s|%s|%s> = %4.9f'%(LS_root[i].get('braStateNumber'),LS_root[i].get('name'),LS_root[i].get('ketStateNumber'),float(LS_root[i].get('value'))))

    else :
        for i in range(len(LS_root)) :
            #if LS_root[i].get('method') == 'MRCI LS_I-I(TOT)' and LS_root[i].get('name') == lstype :
            #if LS_root[i].get('method') in ['MRCI','MRCI LS_I-I(TOT)']  and LS_root[i].get('name') == lstype :
            if LS_root[i].get('method') == method   and LS_root[i].get('name') == lstype :
                index+=1
                if index >= 82:
                    break
                else :
                    # print(index)
                    if 'braStateNumber' in LS_root[i].keys():
                        if 'imaginary'  in LS_root[i].keys():
                            # print(float(LS_root[i].get('value')) * autocminv )
                            so_mat[int(LS_root[i].get('braStateNumber'))-1,int(LS_root[i].get('ketStateNumber'))-1] = \
                                complex(0,float(LS_root[i].get('value')) * au_to_cminv)
                            # so_mat[int(LS_root[i].get('ketStateNumber'))-1,int(LS_root[i].get('braStateNumber'))-1] = \
                            #     np.conj(complex(0,float(LS_root[i].get('value')) * au_to_cminv))
                            # print(int(LS_root[i].get('braStateNumber'))-1,int(LS_root[i].get('ketStateNumber'))-1,float(LS_root[i].get('value')) * au_to_cminv)
                        else:
                            so_mat[int(LS_root[i].get('braStateNumber'))-1,int(LS_root[i].get('ketStateNumber'))-1] = \
                                complex(float(LS_root[i].get('value')) * au_to_cminv,0.0)
                            # so_mat[int(LS_root[i].get('ketStateNumber'))-1,int(LS_root[i].get('braStateNumber'))-1] = \
                            #     complex(-float(LS_root[i].get('value')) * au_to_cminv,0.0)
                        # print(so_mat[int(LS_root[i].get('braStateNumber'))-1,int(LS_root[i].get('ketStateNumber'))-1],so_mat[int(LS_root[i].get('ketStateNumber')) - 1, int(LS_root[i].get('braStateNumber')) - 1])
                    else :
                        if 'imaginary'  in LS_root[i].keys():
                            # print(float(LS_root[i].get('value')) * autocminv)
                            so_mat[int(LS_root[i].get('stateNumber'))-1,int(LS_root[i].get('stateNumber'))-1] = \
                                complex(0,float(LS_root[i].get('value')) * au_to_cminv)
                            # print(int(LS_root[i].get('stateNumber'))-1,int(LS_root[i].get('stateNumber'))-1,float(LS_root[i].get('value')) * au_to_cminv)
                        else:
                            #print(complex(float(LS_root[i].get('value'))* autocminv,0.0))
                            # val = complex(float(LS_root[i].get('value'))*au_to_cminv,0.0)
                            so_mat[int(LS_root[i].get('stateNumber'))-1,int(LS_root[i].get('stateNumber'))-1] = \
                                complex(-float(LS_root[i].get('value'))*au_to_cminv,0.0)
                        # print(so_mat[int(LS_root[i].get('stateNumber')) - 1, int(LS_root[i].get('stateNumber')) - 1])

    # print(so_mat)
    return so_mat

def get_diff_3(mat1,mat2,mat3):
    print('mat1.shape',mat1.shape)
    print('mat2.shape',mat2.shape)
    print('mat3.shape',mat3.shape)
    for j in range(mat1.shape[0]):
        for i in range(mat1.shape[1]):
            # if LSOP_s[n_sub_mat][i][j].imag != 0:
            print("%2d %2d %8.3f +( %8.3f j) %8.3f +( %8.3f j) %8.3f %8.3f"
                  % (i + 1, j + 1,
                     mat1[i][j].real,
                     mat1[i][j].imag,
                     mat2[i][j].real,
                     mat2[i][j].imag,
                     mat3[i][j].real,
                     mat3[i][j].imag))
                     # mat1[i][j].real - mat2[i][j].real,
                     # mat1[i][j].imag - mat2[i][j].imag))

def get_diff(mat1,mat2,op='-',offset1=0,offset2=0,p=False):
    real_error_sum=0
    imag_error_sum=0
    n=0
    r_temp=[]
    i_temp=[]
    print('mat1.shape',mat1.shape)
    print('mat2.shape',mat2.shape)
    if op == 'none':
        for j in range(mat1.shape[1]):
            for i in range(mat1.shape[0]):
                # if LSOP_s[n_sub_mat][i][j].imag != 0:
                print("%2d %2d %8.3f +( %8.3f j) %8.3f +( %8.3f j)" %
                      (i + 1 + offset1, j + 1 + offset2, mat1[i][j].real, mat1[i][j].imag, mat2[i][j].real, mat2[i][j].imag))
    if op == '-':
        for j in range(mat1.shape[1]):
            for i in range(mat1.shape[0]):
                if p == True :
                # if LSOP_s[n_sub_mat][i][j].imag != 0:
                    print("%2d %2d %8.3f +( %8.3f j) %8.3f +( %8.3f j) %8.3f %8.3f"
                      % (i + 1 + offset1, j + 1 + offset2,
                         mat1[i][j].real,
                         mat1[i][j].imag,
                         mat2[i][j].real,
                         mat2[i][j].imag,
                         mat1[i][j].real - mat2[i][j].real,
                         mat1[i][j].imag - mat2[i][j].imag))
                r_temp.append(abs(mat1[i][j].real - mat2[i][j].real))
                i_temp.append(abs(mat1[i][j].imag - mat2[i][j].imag))
                real_error_sum+=abs((mat1[i][j].real - mat2[i][j].real))
                imag_error_sum+=abs((mat1[i][j].imag - mat2[i][j].imag))
                n+=1
        print('Maximum_error:', max(r_temp), max(i_temp))
        print('Average_error:', real_error_sum/n, imag_error_sum/n)

    if op == '+':
        for j in range(mat1.shape[1]):
            for i in range(mat1.shape[0]):
                # if LSOP_s[n_sub_mat][i][j].imag != 0:
                print("%2d %2d %8.3f +( %8.3f j) %8.3f +( %8.3f j) %8.3f %8.3f"
                      % (i + 1 + offset1, j + 1 + offset2,
                         mat1[i][j].real,
                         mat1[i][j].imag,
                         mat2[i][j].real,
                         mat2[i][j].imag,
                         mat1[i][j].real + mat2[i][j].real,
                         mat1[i][j].imag + mat2[i][j].imag))

def get_fac(mat1,mat2,p=False):

    """

    :param mat1:
    :param mat2:
    :return:
    """

    #fac1=mat1[0][1].real / mat2[0][1].real
    # print('Finding factors for the matrix elements')
    if p != True :
        # print('p is False')

        if abs(mat2[0][1].imag) > 1e-8 :
            fac1=mat1[2][5].imag / mat2[2][5].imag
            fac2=mat1[0][1].imag / mat2[0][1].imag
            if fac1 - fac2 > 1e-2 :
                print(fac1-fac2,'is > 1e-3')
                print('something is wrong! ')
            else :
                fac = fac2
        else :
            #print('%8.4f Enconterend divide by 0'%(abs(mat2[0][1].imag)))
            fac = 'inf'

        #if fac1 - fac2 > 0.0001:
        #    print('something is wrong! ')
        #else:
        #    fac = fac2
        #    return fac
        return fac
        print('fac:', fac)
    else :
        print('p is True')
        print("mat1 = %8.4f mat2 = %8.4f"%(mat1[0][1].imag,mat2[0][1].imag))
        for j in range(len(mat1)):
            for i in range(len(mat1)):
                # if LSOP_s[n_sub_mat][i][j].imag != 0:
                print("%2d %2d %8.3f +( %8.3f j) %8.3f +( %8.3f j) %9.6f %9.6f"
                      % (i + 1, j + 1,
                         mat1[i][j].real,
                         mat1[i][j].imag,
                         mat2[i][j].real,
                         mat2[i][j].imag,
                         mat1[i][j].real / mat2[i][j].real,
                         mat1[i][j].imag / mat2[i][j].imag))


def print_mat(mat):
    print('mat.shape',mat.shape)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            print("%d %d %8.4f %8.4f" %(j, i, mat[j][i].real, mat[j][i].imag))

def build_s_sp1_lsmat(LS_org='',S1=3,S2=2,nstates=9,files_path='Example/Na2Co2TeO6',x_prefix='',y_prefix='',z_prefix=''):
    from read_molpro import read_lsop_socE as r_lsop_E

    # nstates = 9
    # S1 = 3
    # S2 = 2
    M1 = 2 * S1 + 1
    M2 = 2 * S2 + 1

    # z_prefix=str(2*S)+str(2*S)+str(2*S)+str(2*S)+'.xml'
    # z_prefix = str(2 * S1) + str(2 * S2) + str(2 * S2) + str(2 * S2) + '_LSZ.xml'
    # x_prefix = str(2 * S1) + str(2 * S2) + str(2 * S1) + str(2 * S2) + '_LSX.xml'
    # y_prefix = str(2 * S1) + str(2 * S2) + str(2 * S1) + str(2 * S2) + '_LSY.xml'
    #
    file = files_path + '/LSX/' + x_prefix
    # print(file)
    LSX = np.array(LS_tran_mat_elem(file, 'LSX', S1, S2))
    # print('LSX.shape', LSX.shape)
    # print_mat(LSX)
    file = files_path + '/LSy/' + y_prefix
    LSY = np.array(LS_tran_mat_elem(file, 'LSY', S1, S2))
    # print('LSY.shape', LSY.shape)

    file = files_path + '/LSZ/' + z_prefix
    # print(file)
    LSZ = np.array(LS_tran_mat_elem(file, 'LSZ', S1, S2))
    # print('LSZ.shape', LSZ.shape)

    LSXY = LSX + LSY
    # print('LSXY.shape', LSXY.shape)

    f = files_path + '/cas10_soc.xml'

    # ls32_split = np.array(blockshaped(ls32, 9, 9))

    new_LS = np.zeros((M1 * nstates, M2 * nstates), dtype=np.complex128)
    print(new_LS.shape)
    new_LS_split = np.array(blockshaped(new_LS, nstates, nstates))
    #
    b = []
    for i in range(0, M2, 1):
        fac = get_fac(LS_org[0], LS_org[i * (M2 + 1)])
        b.append(fac)
        new_LS_split[i * (M2 + 1)] = LSXY / b[i]
    # print(b)

    a = []
    for i in range(0, M2, 1):
        fac = get_fac(LS_org[M2], LS_org[M2 + i * (M2 + 1)])
        a.append(fac)
        new_LS_split[M2 + i * (M2 + 1)] = LSZ / a[i]
    # print(a)

    bs = []
    for i in range(0, M2, 1):
        fac = get_fac(LS_org[0], np.conjugate(LS_org[2 * M2 + i * (M2 + 1)]))
        bs.append(fac)
        new_LS_split[2 * M2 + i * (M2 + 1)] = np.conjugate(LSXY) / bs[i]
    # print(bs)
    # get_fac(ls32_split[0],ls32_split[34],p=True)

    for i in range(M2):
        # print(new_LS[i*9:i*9+9, i*9:i*9+9].shape,new_LS_split[i*(M2+1)].shape)
        new_LS[i * 9:i * 9 + 9, i * 9:i * 9 + 9] = new_LS_split[i * (M2 + 1)]
        new_LS[(i + 1) * 9:(i + 1) * 9 + 9, i * 9:i * 9 + 9] = new_LS_split[M2 + i * (M2 + 1)]
        new_LS[(i + 2) * 9:(i + 2) * 9 + 9, i * 9:i * 9 + 9] = new_LS_split[2 * M2 + i * (M2 + 1)]

    new_ls32_ct = np.transpose(np.conjugate(new_LS))

    LSOP_s_sp1 = new_LS
    LSOP_s_sp1_ct = new_ls32_ct

    return LSOP_s_sp1, LSOP_s_sp1_ct


if __name__ == '__main__':
    import numpy as np
    file = 'Example/Na2Co2TeO6/LSZ/6444_LSZ.xml'
    LSZ = np.array(LS_tran_mat_elem(file, 'LSZ', 3, 2))
    print_mat(LSZ)
    # print(LSZ)

def random():
    import numpy as np
    from read_molpro import split
    from read_molpro import read_mrci_energies as rme



    S=1
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


    for k in range(9):
        l = k + 9
        LSOP_ex[k:l:1, k:l:1] = LSZ + cf_h

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

for i in range(len(mat)):
    for j in range(len(mat)):
        #print("%2d %2d %8.4f %8.4f " %(i,j,mat[i][j].real,mat[i][j].imag,))
        print("%2d %2d %8.4f %8.4f %2d %2d %8.4f %8.4f"
              %(i,j,mat[i][j].real,mat[i][j].imag,i,j,mat2[i][j].real,mat2[i][j].imag))
        # print("%8.4f %8.4f %8.4f" %(LSZ_org[i][j].imag,LSOP_org_1010[i][j].imag,LSZ_org[i][j].imag-LSOP_org_1010[i][j].imag))

mat3=LSOP_org_split[8]

"""





