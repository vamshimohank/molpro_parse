


def get_w3j_fac(filename):

    from read_molpro import read_mrci_energies as rme
    from read_molpro import read_lsop_socE, LS_matrix_elements_xml, split,lsop_read_mod
    import numpy as np


    #au=27.2114
    cminv=219474.63068
    au_to_cminv= cminv

    molpro_file='Example/all_r9.xml'
    m = rme('Example/all_r9.xml')

    #molpro_file='Example/m6_r9.xml'
    #m = rme('Example/m6_r9.xml')



    diag_e = np.array(m['reference Energies'])

    cf_h = np.diag(diag_e-np.min(diag_e)) * au_to_cminv
    #print(cf_h)

    #LSOP,soc_e = read_lsop_socE(molpro_file)
    LSOP,SOC_E = lsop_read_mod(molpro_file)
    print(LSOP.shape)
    #m6lsop1=LSOP[0: 9:1, 0: 9:1]
    m6lsop=LSOP[0: 63:1, 0: 63:1]
    m4lsop=LSOP[63:108:1, 63:108:1]
    m2lsop=LSOP[108:135:1, 108:135:1]
    m0lsop=LSOP[135:144:1, 135:144:1]