


def get_w3j_fac(filename):

    from read_molpro import read_mrci_energies as rme
    from read_molpro import LS_matrix_elements_xml, split,lsop_read_mod
    import numpy as np


    au_to_cminv = 219474.63068

    molpro_file=filename
    # m = rme(molpro_file)

    LSOP,SOC_E = lsop_read_mod(molpro_file)

    m6lsop = LSOP[0: 63:1, 0: 63:1]     # 2S=6 part of the SOC matrix shape(63,63)
    m4lsop = LSOP[63:108:1, 63:108:1]   # 2S=4 part of the SOC matrix shape(45,45)
    m2lsop = LSOP[108:135:1, 108:135:1] # 2S=2 part of the SOC matrix shape(27,27)
    m0lsop = LSOP[135:144:1, 135:144:1] # 2S=0 part of the SOC matrix shape(9,9)

    # split each Spin Multiplet into blocks of 9 states
    m6_LSOP_r9 = np.array(split(m6lsop, 9, 9))
    m4_LSOP_r9 = np.array(split(m4lsop, 9, 9))
    m2_LSOP_r9 = np.array(split(m2lsop, 9, 9))
    m0_LSOP_r9 = np.array(split(m0lsop, 9, 9))

    m6_d_ms_1 = m6_LSOP_r9[0][2][1].imag/m6_LSOP_r9[8][2][1]
    m6_d_ms_2 = m6_LSOP_r9[0][2][1].imag/m6_LSOP_r9[16][2][1]
    print(m6_d_ms_1,m6_d_ms_2)

    m4_d_ms_1 = m4_LSOP_r9[0][2][1].imag/m4_LSOP_r9[6][2][1]
    m4_d_ms_2 = m4_LSOP_r9[0][2][1].imag/m4_LSOP_r9[12][2][1]
    print(m4_d_ms_1,m4_d_ms_2)

    m2_d_ms_1 = m2_LSOP_r9[0][2][1].imag/m2_LSOP_r9[4][2][1]
    print(m2_d_ms_1)

def build_LSOP_from_LS_mat_elem(molpro_file):

    from read_molpro import LS_matrix_elements_xml

    LSZ=LS_matrix_elements_xml(molpro_file,'LSZ')

    print(LSZ.shape)
    return LSZ
    # for i in range(len(LSZ)):
    #     for j in range(len(LSZ)):
    #         print("%d %d %5.3f" % (i,j,LSZ[i][j]))
    #LSOP_1_1_ex=LSZ+cf_h

    LSX=LS_matrix_elements_xml(molpro_file,'LSX')
# LSY=LS_matrix_elements_xml(molpro_file,'LSY')

#get_w3j_fac('all_r9.xml')
LSZ=build_LSOP_from_LS_mat_elem('all_r9.xml')
for i in range(len(LSZ)):
        for j in range(len(LSZ)):
            print("%d %d %8.3f %8.3f" % (i,j,LSZ[j][j].real,LSZ[i,j].imag))