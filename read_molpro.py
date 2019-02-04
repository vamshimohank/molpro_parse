#!/usr/bin/python


#Pilot implementation of XAS and RIXS intensities calculation for Li2CuO2

# The LSOP and WFs extraction parts initially done by Vamshi M. Katukuri (2013)
# Improved and expanded by Nikolay A. Bogdanov (2014)
#
# The RIXS and XAS parts by Nikolay A. Bogdanov (2014, 2015)
# Licensed under the GNU General Public License v3.0


#filename = sys.argv[1]
#typ = sys.argv[2]

#filename = 'tm.out'
#typ='rixs_cas'

#f = open(filename).readlines()                     # read Molpro output file
#tst = open('abs_cross.txt','w')

def split(array, nrows, ncols):
    """Split a matrix into sub-matrices."""

    r, h = array.shape
    return (array.reshape(h//nrows, nrows, -1, ncols)
                 .swapaxes(1, 2)
                 .reshape(-1, nrows, ncols))


def read_wfc(f,write_files=False):
    import numpy as np
    wfc_lookup = 'Eigenvectors of spin-orbit matrix'
    wfc_lookup_end = ' Composition of spin-orbit eigenvectors'
    WFstart_line=[]

    for num1, line in enumerate(f):                                     # Find all the blocks we want to extract
        if wfc_lookup in line:
            WFstart_line.append(num1+7)                                 # Start of the wave function coefficients
    WF=[]
    nstates=[]
    for i, num in enumerate(WFstart_line):
        nstates.append(int(f[num - 10][0:4]))  # last line of energies output
        WF.append(np.zeros((nstates[i], nstates[i]), dtype=np.complex128))  # fill WFs array with zeros
        WF[i] = WF[i] + 100 + 100j  # fill WFs array with 100's to prevent errors

    #
    #print (nstates)



    for ns, n_states in enumerate(nstates):
        if write_files:
            wfc_fname = typ + ".wfc." + str(ns) + "." + str(n_states) + "s"  # create files for WFs
            wfc_out = open(wfc_fname, 'w')

    #############
    ###### WFC block details
    #############
        wfc_bcs_d = 8  # way it is printed in molpro
        wfc_block_col_size = 8  # equal to default if n_states < default
        if n_states <= wfc_block_col_size:                              # if we have less than 8 SO-states
            wfc_n_blocks = 1
            wfc_block_col_size = n_states
            wfc_last_block_col_size = 0
        else:
            wfc_n_blocks= int(n_states/wfc_block_col_size)                   # find out how many 8-column blocks are printed
            if wfc_n_blocks*wfc_block_col_size != n_states:             # if number of SO-states isn't 8*N
                wfc_last_block_col_size = n_states-(wfc_n_blocks*wfc_block_col_size) # add last block with width <8
                wfc_n_blocks=wfc_n_blocks+1
            else:
                wfc_block_col_size = 8                                  # if number of SO-states is 8*N
                wfc_last_block_col_size = 0                             # no additional blocks needed
        wfc_block_rows = n_states*3                            # row size for each block
        wfc_block_gap = 3

    #############
    ###### Extraction of WFC coefficients
    #############

        for nbl in range(0,wfc_n_blocks):                               # loop over the No of blocks
            if wfc_last_block_col_size != 0 and nbl == wfc_n_blocks-1:  # make less columns for last block
                wfc_block_col_size = wfc_last_block_col_size            #
    #                                                                   # specify block of 2*Nroots rows and 8 columns
            for n in range(WFstart_line[ns] + nbl*wfc_block_rows + nbl*wfc_block_gap, WFstart_line[ns] + (nbl+1)*wfc_block_rows + nbl*wfc_block_gap) :
                    if f[n]!='\n' and f[n]!=' \n' :                     # not to consider empty lines
                        if f[n][2]!=' ' :                               # rows with real parts also contain No, and (j, mj) information (only No is now used)
                            for col in range(0,wfc_block_col_size) :
                                                                        # select only relevant part of the line and read
                                WF[ns][int(f[n].split()[0])-1,col+nbl*wfc_bcs_d]=complex(float(f[n][22:135].split()[col]),0)
                        else :                                          # rows with imaginary parts have spaces in the beginning
                            for col in range(0,wfc_block_col_size) :
    #                                                                   #  col+nbl*8 is 'absolute' No of column
                                WF[ns][int(f[n-1].split()[0])-1][col+nbl*wfc_bcs_d]=WF[ns][int(f[n-1].split()[0])-1][col+nbl*wfc_bcs_d]+complex(0,float(f[n][22:135].split()[col]))
        for nr in range(n_states) :                                     # print to file
            prnt=''                                                     # specify empty line, then 'collect' the whole row
            for nc in range(n_states):
                prnt += str("%11.9f"%WF[ns].real[nr,nc])+str("+I(%11.9f) "%WF[ns].imag[nr,nc]) # normal output, WFs in columns
    #            prnt += str("%11.9f"%WF[ns].real[nc,nr])+str("+I(%11.9f) "%WF[ns].imag[nc,nr]) # transposed output, WFs in rows
            #print >> wfc_out, prnt
            if write_files :
                print(prnt,file=wfc_out)
        if write_files :
            wfc_out.close()

    return nstates,WF

def lsop_read_mod(filename,write_files=False):
    import numpy as np

    f = open(filename).readlines()


    for num1, line in enumerate(f):  # Find all the blocks we want to extract
        if 'VMK Spin-Orbit Matrix (CM-1)' in line:
            LSOPstart_line=num1 + 3  # Start of the Spin-Orbit Matrix data
            LSOPlines = True
        #        print('Found LSOP at line', num1)
            print(LSOPstart_line)

        if 'No symmetry adaption' in line and LSOPlines == True:
            LSOPend_line = num1 - 1  # End of the Spin-Orbit Matrix data
            print('LSOPend_line=',LSOPend_line)
        if 'Spin-orbit eigenstates' in line and LSOPlines == True:
            if 'No symmetry adaption' in f[num1-3] :
                continue
            else:
                LSOPend_line = num1 - 3  # End of the Spin-Orbit Matrix data
                print('LSOPend_line=',LSOPend_line)

        if 'Nr         E' in line :
            soc_e_start_line=num1 + 2

        if 'Eigenvectors of spin-orbit matrix' in line:
            if '=====' in f[num1+1]:
                WFstart_line=num1 + 7
                print(WFstart_line)

    nstates=int(f[WFstart_line - 10][0:4])
    #print(nstates)
    SOC_E = np.zeros(nstates,dtype=float)

    for n in range(soc_e_start_line, WFstart_line - 9):
        SOC_E[int(f[n].split()[0]) - 1] = float(f[n].split()[1])  # *hartree
        if write_files:
            print(SOC_E[ns][int(f[n].split()[0]) - 1], file=socE_out)
    bra=[]
    ket=[]
    braS=[]
    braMs=[]
    ket={}
    LSOP=np.zeros((nstates,nstates),dtype=np.complex128)
    for num1, line in enumerate(f):  # Find all the blocks we want to extract
        if num1 >= LSOPstart_line and num1 <= LSOPend_line and line != '\n' :
            # print(line.split())
            bra=int(line.split()[0])
            ket=int(line.split()[10])
            LSOP[bra-1][ket-1]=np.complex(float(line.split()[11]),float(line.split()[13]))

    return LSOP,SOC_E


def read_lsop_socE_org(filename,write_files=False):
    import numpy as np

    f = open(filename).readlines()

    WFstart_line=[]
    LSOPstart_line = []
    LSOPend_line = []
    LSOPlines = False


    for num1, line in enumerate(f):  # Find all the blocks we want to extract

        if 'Eigenvectors of spin-orbit matrix' in line:
            WFstart_line.append(num1 + 7)

        if ' Spin-Orbit Matrix (CM-1)' in line:
            LSOPstart_line.append(num1 + 5)  # Start of the Spin-Orbit Matrix data
            LSOPlines = True
        #        print('Found LSOP at line', num1)
        # if 'No symmetry adaption' in line and LSOPlines == True:
        #     LSOPend_line.append(num1 - 2)  # End of the Spin-Orbit Matrix data
        #     LSOPlines = False
        if 'Spin-orbit eigenstates' in line and LSOPlines == True:
            LSOPend_line.append(num1 - 3)  # End of the Spin-Orbit Matrix data
            LSOPlines = False
            print(LSOPend_line)
    # if LSOPlines==True or WFlines==True or  NOCIlines==True :
    #if LSOPlines == True or NOCIlines == True:
    #     print('One of the blocks has start, but has no end')
    #
    if len(LSOPstart_line) != len(WFstart_line):
        print('# WF blocks found is', len(WFstart_line), 'and of # L.S operator blocks is', len(LSOPstart_line),
              'Something is wrong! check the output file...')
    nstates = []
    LSOP = []
    SOC_E = []
    for i, num in enumerate(WFstart_line):
        nstates.append(int(f[num - 10][0:4]))
        LSOP.append(np.zeros((nstates[i], nstates[i]), dtype=np.complex128))  #
        LSOP[i] = LSOP[i] + 100 + 100j  #
        SOC_E.append(np.zeros(nstates[i]))  # fill energies arrays with zeros

        #so_mat = np.matlib.zeros((nstates[i], nstates[i]), dtype=np.complex128)

    for ns, n_states in enumerate(nstates):
        if write_files:
            lsop_fname = typ + ".lsop." + str(ns) + "." + str(n_states) + "s"  # create files for L.S operators
            lsop_out = open(lsop_fname, 'w')
            socE_fname = typ + ".esoc." + str(ns) + "." + str(n_states) + "s"  # create files for energies
            socE_out = open(socE_fname, 'w')

        #############
        ###### LSOP block details
        #############
        lop_block_rows = n_states * 3  # row size for each block
        lop_block_gap = 3
        lsop_bcs_d = 10  # way it is printed in molpro
        lsop_block_col_size = lsop_bcs_d  # equal to default if n_states < default
        if n_states <= lsop_block_col_size:
            lsop_n_blocks = 1
            lsop_block_col_size = n_states
            lsop_last_block_col_size = 0
        else:
            lsop_n_blocks = int(n_states / lsop_block_col_size)  # No of blocks the matrix is written
            if lsop_n_blocks * lsop_block_col_size != n_states:
                lsop_last_block_col_size = n_states - (lsop_n_blocks * lsop_block_col_size)
                lsop_n_blocks = lsop_n_blocks + 1
            else:
                lsop_last_block_col_size = 0

        #############
        ###### LSOP extraction
        #############
        #
        for nbl in range(0, lsop_n_blocks):  # loop over the No of blocks
            if lsop_last_block_col_size != 0 and nbl == lsop_n_blocks - 1:
                lsop_block_col_size = lsop_last_block_col_size
                #                                                                   # loop over block lines
            for n in range(LSOPstart_line[ns] + nbl * lop_block_rows + nbl * lop_block_gap,
                           LSOPstart_line[ns] + (nbl + 1) * lop_block_rows + nbl * lop_block_gap):
                if f[n] != '\n' and f[n] != ' \n':  # not to consider empty lines
                    if f[n][
                        2] != ' ':  # rows with real parts also contain No, and (j, mj) information (only No is now used)
                        for col in range(0, lsop_block_col_size):
                            LSOP[ns][int(f[n].split()[0]) - 1][col + nbl * lsop_bcs_d] = complex(
                                float(f[n][19 + 11 * col:30 + +11 * col]),
                                0)  # select only relevant part of the line and read
                    else:  # rows with imaginary parts have spaces in the beginning
                        for col in range(0, lsop_block_col_size):
                            LSOP[ns][int(f[n - 1].split()[0]) - 1][col + nbl * lsop_bcs_d] = \
                            LSOP[ns][int(f[n - 1].split()[0]) - 1][col + nbl * lsop_bcs_d] + complex(0, float(f[n][
                                                                                                              19 + 11 * col:30 + 11 * col]))  # 19 symbol in line is the first for LSOP, every number is 11 symbols => 30 - end of the first number
        #
        if write_files :
            for nr in range(n_states):  # print to file
                prnt = ''  # specify empty line, then 'collect' the whole row
                for nc in range(n_states):
                    prnt += str(LSOP[ns].real[nr, nc]) + "+I(" + str(
                        LSOP[ns][nr][nc].imag) + ") "  # normal output, same way as in molpro
                #            prnt += str("%11f"%LSOP[ns][nr][nc].real)+str("+I(%11f) "%LSOP[ns].imag[nr,nc]) # normal output, same way as in molpro (other format(?))
                #            prnt += str("%11f"%LSOP[ns].real[nc,nr])+str("+I(%11f) "%LSOP[ns].imag[nc,nr]) # transposed output, WFs in rows

                print(prnt, file=lsop_out)
            lsop_out.close()

        #############
        ###### Energies extraction
        #############
        #
        for n in range(LSOPend_line[ns] + 10, WFstart_line[ns] - 9):
            SOC_E[ns][int(f[n].split()[0]) - 1] = float(f[n].split()[1])  # *hartree
            if write_files :
                print(SOC_E[ns][int(f[n].split()[0])-1],file=socE_out)
        if write_files :
            socE_out.close()
    LSOP = np.array(LSOP)
    return LSOP,SOC_E



def read_lsop_socE(filename,write_files=False):
    import numpy as np

    f = open(filename).readlines()

    WFstart_line=[]
    LSOPstart_line = []
    LSOPend_line = []
    LSOPlines = False


    for num1, line in enumerate(f):  # Find all the blocks we want to extract

        if 'Eigenvectors of spin-orbit matrix' in line:
            WFstart_line.append(num1 + 7)

        if ' Spin-Orbit Matrix (CM-1)' in line:
            LSOPstart_line.append(num1 + 5)  # Start of the Spin-Orbit Matrix data
            LSOPlines = True
        #        print('Found LSOP at line', num1)
        # if 'No symmetry adaption' in line and LSOPlines == True:
        #     LSOPend_line.append(num1 - 2)  # End of the Spin-Orbit Matrix data
        #     LSOPlines = False
        #     print(LSOPend_line)
        if 'Spin-orbit eigenstates' in line and LSOPlines == True:
            LSOPend_line.append(num1 - 4)  # End of the Spin-Orbit Matrix data
            LSOPlines = False
            print(LSOPend_line)
    # if LSOPlines==True or WFlines==True or  NOCIlines==True :
    #if LSOPlines == True or NOCIlines == True:
    #     print('One of the blocks has start, but has no end')
    #
    if len(LSOPstart_line) != len(WFstart_line):
        print('# WF blocks found is', len(WFstart_line), 'and of # L.S operator blocks is', len(LSOPstart_line),
              'Something is wrong! check the output file...')
    nstates = []
    LSOP = []
    SOC_E = []
    for i, num in enumerate(WFstart_line):
        nstates.append(int(f[num - 10][0:4]))
        LSOP.append(np.zeros((nstates[i], nstates[i]), dtype=np.complex128))  #
        LSOP[i] = LSOP[i] + 100 + 100j  #
        SOC_E.append(np.zeros(nstates[i]))  # fill energies arrays with zeros

        #so_mat = np.matlib.zeros((nstates[i], nstates[i]), dtype=np.complex128)

    for ns, n_states in enumerate(nstates):
        if write_files:
            lsop_fname = typ + ".lsop." + str(ns) + "." + str(n_states) + "s"  # create files for L.S operators
            lsop_out = open(lsop_fname, 'w')
            socE_fname = typ + ".esoc." + str(ns) + "." + str(n_states) + "s"  # create files for energies
            socE_out = open(socE_fname, 'w')

        #############
        ###### LSOP block details
        #############
        lop_block_rows = n_states * 3  # row size for each block
        lop_block_gap = 3
        lsop_bcs_d = 10  # way it is printed in molpro
        lsop_block_col_size = lsop_bcs_d  # equal to default if n_states < default
        if n_states <= lsop_block_col_size:
            lsop_n_blocks = 1
            lsop_block_col_size = n_states
            lsop_last_block_col_size = 0
        else:
            lsop_n_blocks = int(n_states / lsop_block_col_size)  # No of blocks the matrix is written
            if lsop_n_blocks * lsop_block_col_size != n_states:
                lsop_last_block_col_size = n_states - (lsop_n_blocks * lsop_block_col_size)
                lsop_n_blocks = lsop_n_blocks + 1
            else:
                lsop_last_block_col_size = 0

        #############
        ###### LSOP extraction
        #############
        #
        for nbl in range(0, lsop_n_blocks):  # loop over the No of blocks
            if lsop_last_block_col_size != 0 and nbl == lsop_n_blocks - 1:
                lsop_block_col_size = lsop_last_block_col_size
                #                                                                   # loop over block lines
            for n in range(LSOPstart_line[ns] + nbl * lop_block_rows + nbl * lop_block_gap,
                           LSOPstart_line[ns] + (nbl + 1) * lop_block_rows + nbl * lop_block_gap):
                print(n)
                if f[n] != '\n' and f[n] != ' \n':  # not to consider empty lines
                    if f[n][2] != ' ':  # rows with real parts also contain No, and (j, mj) information (only No is now used)
                        for col in range(0, lsop_block_col_size):
                            LSOP[ns][int(f[n].split()[0]) - 1][col + nbl * lsop_bcs_d] = complex(
                                float(f[n][19 + 11 * col:30 + +11 * col]),
                                0)  # select only relevant part of the line and read
                    else:  # rows with imaginary parts have spaces in the beginning
                        for col in range(0, lsop_block_col_size):
#                            print(LSOP[ns]) #[int(f[n-1])])
                            print(n,f[n - 1].split())
                            LSOP[ns][int(f[n - 1].split()[0]) - 1][col + nbl * lsop_bcs_d] = \
                            LSOP[ns][int(f[n - 1].split()[0]) - 1][col + nbl * lsop_bcs_d] + \
                            complex(0, float(f[n][19 + 11 * col:30 + 11 * col]))  # 19 symbol in line is the first for LSOP, every number is 11 symbols => 30 - end of the first number
        #
        if write_files :
            for nr in range(n_states):  # print to file
                prnt = ''  # specify empty line, then 'collect' the whole row
                for nc in range(n_states):
                    prnt += str(LSOP[ns].real[nr, nc]) + "+I(" + str(
                        LSOP[ns][nr][nc].imag) + ") "  # normal output, same way as in molpro

                print(prnt, file=lsop_out)
            lsop_out.close()

        #############
        ###### Energies extraction
        #############
        #
        for n in range(LSOPend_line[ns] + 10, WFstart_line[ns] - 9):
            SOC_E[ns][int(f[n].split()[0]) - 1] = float(f[n].split()[1])  # *hartree
            if write_files :
                print(SOC_E[ns][int(f[n].split()[0])-1],file=socE_out)
        if write_files :
            socE_out.close()
    LSOP = np.array(LSOP)
    return LSOP,SOC_E
#

def read_noci_tm(f,write_files=False):
    import numpy as np

    NOCIstart_line = []
    NOCIend_line = []
    NOCIlines = False

    NOCI_lookup_start = ' Bi-orthogonal integral transformation finished'
    NOCI_lookup_end = ' **********************************************************************************************************************************'

    for num1, line in enumerate(f):
        if NOCI_lookup_start in line:
            NOCIstart_line.append(num1 + 1)
            NOCIlines = True
        if NOCI_lookup_end in line and NOCIlines == True:
            NOCIend_line.append(num1 - 2)
            NOCIlines = False

    dNCI = []
    D = [[] for ni in range(3)]  # DX=D[0], DY=D[1], DZ=D[2]
    dN = 0

    for nN in range(len(NOCIstart_line)):
        for lin in f[NOCIstart_line[nN]:NOCIend_line[nN]]:
            for char in ['<', '>', '|']:
                lin = lin.replace(char, ' ')
            if lin[0] != '\n' and lin.split()[0] == '!MRCI' and lin.split()[3] == 'H':
                dN = max(dN, float(lin.split()[2]), float(lin.split()[4]))
        dN = int(dN - 0.1)
        dNCI.append(dN)
        for x in range(3):
            D[x].append(np.zeros((dN, dN)))

    #############
    ###### Reading and writing dipole moment matrices from NOCI blocks
    #############
    #
    for nN, dN in enumerate(dNCI):
        if write_files :
            dm_bname = ".dm." + str(nN) + "." + str(dN) + "s."  # base name for dipole moments files
            fileDX = open(typ + dm_bname + "X", 'w')  # dm.0.8s.X_ci
            fileDY = open(typ + dm_bname + "Y", 'w')
            fileDZ = open(typ + dm_bname + "Z", 'w')
        #
        for lin in f[NOCIstart_line[nN]:NOCIend_line[nN]]:
            for char in ['<', '>', '|']:  # replace brackets with space
                lin = lin.replace(char, ' ')
            if lin[0] != '\n' and lin.split()[0] == '!MRCI' and lin.split()[3][0:2] == 'DM':  # find DM lines
                if lin.split()[3] == 'DMX':
                    D[0][nN][int(float(lin.split()[2])) - 1][int(float(lin.split()[4])) - 1] = lin.split()[
                        5]  # save matrix element
                elif lin.split()[3] == 'DMY':
                    D[1][nN][int(float(lin.split()[2])) - 1][int(float(lin.split()[4])) - 1] = lin.split()[5]
                elif lin.split()[3] == 'DMZ':
                    D[2][nN][int(float(lin.split()[2])) - 1][int(float(lin.split()[4])) - 1] = lin.split()[5]
                else:
                    print('Something is wrong !MRCI <i.1|DM_|j.1> lines, at line:', lin)
        if write_files :
            for item in D[0][nN]:
                print(item[0], " ".join(map(str, item[1:])),file=fileDX)
            for item in D[1][nN]:
                print(item[0], " ".join(map(str, item[1:])),file=fileDY)
            for item in D[2][nN]:
                print(item[0], " ".join(map(str, item[1:])),file=fileDZ)
    return D

def read_mrci_energies(file_name):
    """

    :param file_name: name of the molpro output file
    :return: multiplicity (int), nstates (int), ref_energies (numpy array), mrci_energies (numpy array), mrci_energies_davidson_fixed (numpy array), mrci_eneriges_davidson_relaxed (numpy array)
    """
    import numpy as np
    f=open(file_name, 'r')
    for line in f:
        if 'Multireference internally contracted CI' in line :
            for lin  in f:
                if 'Number of optimized states:' in lin:
                    nstates = int(lin.split()[4])
                if 'Reference symmetry:' in lin:
                    multiplicity = int(lin.split()[2])
                if 'State     Reference Energy' in lin:
                    ref_energies = []
                    mrci_energies = []
                    mrci_energies_davidson_fixed = []
                    mrci_energies_davidson_relax = []
                    for i, li in enumerate(f):

                        if li[0] != '\n':
                            if i <= nstates :
                                ref_energies.append(float(li.split()[1]))
                            else :
                                ref_energies=np.array(ref_energies)

                                if '!MRCI STATE ' in li and 'Energy' in li :
                                    mrci_energies.append(float(li.split()[4]))
                                if 'Davidson, fixed reference' in li :
                                    mrci_energies_davidson_fixed.append(float(li.split()[3]))
                                if 'Davidson, relaxed reference' in li :
                                    mrci_energies_davidson_relax.append(float(li.split()[3]))

                    mrci_energies=np.array(mrci_energies)
                    mrci_energies_davidson_fixed = np.array(mrci_energies_davidson_fixed)
                    mrci_energies_davidson_relax = np.array(mrci_energies_davidson_relax)

                    return {'multiplicity' : multiplicity, 'nstates' : nstates,
                            'reference Energies' : ref_energies,
                            'MRCI Energies' : mrci_energies,
                            'MRCI Energies with Davidson Fixed frame' : mrci_energies_davidson_fixed,
                            'MRCI Energies with Davidson Relax frame' : mrci_energies_davidson_relax
                            }
                    break

#def read_mrci_energies_xml(file_name):




def LS_matrix_elements_xml(file,lstype,S):


    au=27.2114
    cminv=8.07
    #autocminv= au * cminv * 1000
    autocminv=219474.6308408

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

    # so_mat=numpy.matlib.zeros((nstates,nstates),dtype=np.complex128)
    so_mat=numpy.zeros((nstates,nstates),dtype=np.complex128)
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

    #print(LS_root[0].attrib)
    # print(len(LS_root))

    if 'ECP' in lstype :
        print('In ECP')
        for i in range(len(LS_root)):
            if LS_root[i].get('name') == lstype and LS_root[i].get('method') == 'MRCI LS_I-I(ECP)':

                print('<%s|%s|%s> = %4.9f'%(LS_root[i].get('braStateNumber'),LS_root[i].get('name'),LS_root[i].get('ketStateNumber'),float(LS_root[i].get('value'))))

    else :
        for i in range(len(LS_root)):
            #if LS_root[i].get('method') == 'MRCI LS_I-I(TOT)' and LS_root[i].get('name') == lstype :
            if LS_root[i].get('method') in ['MRCI','MRCI LS_I-I(TOT)']  and LS_root[i].get('name') == lstype :
                if 'braStateNumber' in LS_root[i].keys():
                    if 'imaginary'  in LS_root[i].keys():
                        # print(float(LS_root[i].get('value')) * autocminv )
                        so_mat[int(LS_root[i].get('braStateNumber'))-1,int(LS_root[i].get('ketStateNumber'))-1] = \
                            complex(0,float(LS_root[i].get('value')) * autocminv)
                        so_mat[int(LS_root[i].get('ketStateNumber'))-1,int(LS_root[i].get('braStateNumber'))-1] = \
                            np.conj(complex(0,float(LS_root[i].get('value')) * autocminv))
                    else:
                        so_mat[int(LS_root[i].get('braStateNumber'))-1,int(LS_root[i].get('ketStateNumber'))-1] = \
                            complex(float(LS_root[i].get('value')) * autocminv,0.0)
                        so_mat[int(LS_root[i].get('ketStateNumber'))-1,int(LS_root[i].get('braStateNumber'))-1] = \
                            complex(-float(LS_root[i].get('value')) * autocminv,0.0)
                else :
                    if 'imaginary'  in LS_root[i].keys():
                        # print(float(LS_root[i].get('value')) * autocminv)
                        so_mat[int(LS_root[i].get('stateNumber'))-1,int(LS_root[i].get('stateNumber'))-1] = \
                            complex(0,float(LS_root[i].get('value')) * autocminv)
                    else:
                        #print(complex(float(LS_root[i].get('value'))* autocminv,0.0))
                        val = complex(float(LS_root[i].get('value'))*autocminv,0.0)
                        so_mat[int(LS_root[i].get('stateNumber'))-1,int(LS_root[i].get('stateNumber'))-1] = \
                            complex(-float(LS_root[i].get('value'))*autocminv,0.0)

    #print(so_mat)
    return so_mat

if __name__ == '__main__':
# file='SOC_XML/tranlsx.xml'
# LS_matrix_elements(file,'ECPLSX')
# file='SOC_XML/tranlsy.xml'
# LS_matrix_elements(file,'ECPLSY')
# file='SOC_XML/tranlsz.xml'
# LS_matrix_elements(file,'ECPLSZ')
#file='SOC_XML/lsop_mat.xml'

# LSX=LS_matrix_elements('SOC_XML/lsop_mat_lsx.xml','LSX')
# LSY=LS_matrix_elements('SOC_XML/lsop_mat_lsy.xml','LSY')
# LSZ=LS_matrix_elements('SOC_XML/lsop_mat_lsz.xml','LSZ')

    import numpy as np
    LSOP,soc_e = lsop_read_mod('Example/all_r9.xml')

#print(type(LSOP))
    LSOP9x9=LSOP[0:9:1,0:9:1]
#print(LSOP9x9)
    S=3.0
    LSZ=LS_matrix_elements_xml('Example/all_r9.xml','LSZ',S)
#print(LSZ)
# for i in range(9):
#     for j in range(9):
#         print('<{:2d}|LSZ|{:2d}> = {:.2f}'.format(j+1,i+1,LSZ[j,i]),
#               ' <{:2d}|LSZ|{:2d}> = {:.2f}'.format(j+1,i+1,LSOP9x9[j, i]))
    for i in range(9):
        for j in range(9):
            # print(i,j,'{:.2f} {:.2f}'.format(LSOP9x9[i, j],np.transpose(LSZ[i, j])))
            if LSOP9x9[i,j].real-LSZ[i,j].real >0.4 or LSOP9x9[i,j].imag-LSZ[i,j].imag > 0.4:

                print(i,j,'{:0.2f} {:0.2f}'.
                      format(LSOP9x9[i,j].real-LSZ[i,j].real,
                             LSOP9x9[i,j].imag-LSZ[i,j].imag))


# H_LS=LSX+LSY+LSZ

# print(H_LS)


# file_name='/Users/katukuri/Downloads/r1_9.out'
# multiplicity, nstates, ref_energies, mrci_energies, mrci_energies_davidson_fixed, mrci_energies_davidson_relax = read_mrci_energies(file_name)
#
# print("Multiplicity: ",multiplicity)
# mrci_rel_energies = mrci_energies-min(mrci_energies)
# mrci_rel_energies_davidson_fixed = mrci_energies_davidson_fixed - min(mrci_energies_davidson_fixed)
# mrci_rel_energies_davidson_relax = mrci_energies_davidson_relax - min(mrci_energies_davidson_relax)
# ref_rel_energies=ref_energies-min(ref_energies)
# for i in range(nstates):
#     print('%4.8f  %4.4f  %4.8f  %4.4f  %4.4f  %4.4f  %4.4f  %4.4f'%(ref_energies[i], ref_rel_energies[i]*27.2114,
#                                                                     mrci_energies[i], mrci_rel_energies [i]*27.2114,
#                                                                     mrci_energies_davidson_fixed[i], mrci_rel_energies_davidson_fixed[i]*27.2114,
#                                                                     mrci_energies_davidson_relax[i], mrci_rel_energies_davidson_relax[i]*27.2114))

#D=read_noci_tm(f)
#print(D)
# wf=read_wfc(f)
# print(wf)
# lsop,soce=read_lsop_socE(f)
# for n in range(len(soce)):
#     print(soce[n])


