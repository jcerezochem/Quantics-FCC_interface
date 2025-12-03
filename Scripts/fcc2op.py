#!/usr/bin/env python
import numpy as np
import sys

# Constants and conversions
evtown = 8065.5446811132
autoev = 27.2113961
autown = 2.1947463068e5

def get_state_params(fname):
    f = open(fname)

    model='none'
    for line in f:
        if line[:14] == 'MODEL        =':
            model = line.split('=')[1].strip()
            break

    # Move pin to State 1
    for line in f:
        if '==> State 1 <==' in line:
            break

    # Read freqs
    for line in f:
        if 'FREQUENCIES (cm-1)' in line:
            break
    f.readline() # skip header line
    line = f.readline() # First item
    omega = []
    while not '----' in line:
        w = float(line.split()[1]) / autown
        omega.append(w)
        line = f.readline()
    omega = np.array(omega)
    nvib = len(omega)

    # Move pin to State 2
    for line in f:
        if '==> State 2 <==' in line:
            break

    if model == 'VG' or model == 'VH':
        # Read freqs
        for line in f:
            if 'FREQUENCIES (cm-1)' in line:
                break
        f.readline() # skip header line
        line = f.readline() # First item
        omega_ES = []
        while not '----' in line:
            w = float(line.split()[1]) / autown
            omega_ES.append(w)
            line = f.readline()
        omega_ES = np.array(omega_ES)
        nvib_ = len(omega_ES)
        if nvib != nvib_:
            print('ERROR: Inconsitent vibrational spaces')
            sys.exit()
        # Read gradient
        for line in f:
            if 'State2 Gradient (Q1' in line and '-space) x1e5' in line:
                break
        f.readline() # skip header line
        line = f.readline() # First item
        gradQ = []
        while not '----' in line:
            g = float(line.split()[1]) * 1.e-5
            gradQ.append(g)
            line = f.readline()
        gradQ = np.array(gradQ)
        gradq = gradQ / np.sqrt(omega)
    else:
        print('ERROR: only VG and VH models supported')
        sys.exit()

    # Now get transition properties
    ## Find location of Duschinsky file
    dusch_file = None
    for line in f:
        if 'Printing Duschinski matrix to' in line:
            dusch_file = line.strip().split()[-1].replace("'","")
            break
    ## Transition energy (in eV)
    for line in f:
        if 'Extrapolated vertical energy' in line:
            DE = float(line.split()[3])
            break
    # ZPE (in eV)
    for line in f:
        if ' E01= ' in line:
            ZPE1,ZPE2 = [float(line.split()[i])/evtown for i in [1,3]]

    f.close

    # If VH model, read Dusch file
    if model == 'VH':
        try:
            J = np.loadtxt(dusch_file).reshape([nvib,nvib])
        except:
            print('ERROR: cannot read {dusch_file}')
            sys.exit()

        # Get Hq
        HQ = J @ np.diag(omega_ES**2) @ J.transpose()
        Hq = np.zeros([nvib,nvib])
        for i in range(nvib):
            Hq[i,i] = HQ[i,i] / omega[i]
            for j in range(i):
                Hq[i,j] = HQ[i,j] / np.sqrt(omega[i]*omega[j])
                Hq[j,i] = Hq[i,j]
        Hq -= np.diag(omega)
    else:
        Hq = None

    return omega, DE-ZPE1, gradq, Hq

def write_op(fop, omega, Ener, Gradq, Hqs, zero=1.e-5):
    # Write operator file
    fop = open(fname_out,'w')
    header='''OP_DEFINE-SECTION
 title
 FCclasses VG Model
 end-title
end-op_define-section

PARAMETER-SECTION

# frequencies'''
    print(header, file=fop)

    # Print freqs in eV
    nvib = len(omega)
    omega *= autoev
    for i,w in enumerate(omega):
        print(f'  w{i+1:03g}  =     {w:12.8f}  ,  ev', file=fop)

    # Add energy of the state, which should be DE - ZPE (Quantics sets the Zero at ZPE)
    section='\n# vertical energies (from ZPE) and off diag coupling const (if present)\n'
    for state,e in enumerate(Ener,1):
        section += f'  delta{state:02g}_{state:02g} =        {e}  ,   ev\n'

    print(section, file=fop)

    # Print grad in eV
    print('# linear intrastate parameters', file=fop)
    for state,gradq in enumerate(Gradq):
        gradq *= autoev
        for i,g in enumerate(gradq):
            if np.abs(g) > zero:
                print(f'  KD{state:02g}_{state:02g}_{i+1:03g}  =     {g:12.8f}  ,  ev', file=fop)

    # If VH, then modify ES forces
    for state,Hq in enumerate(Hqs,1):
        if Hq is not None:
            Hq *= autoev / 2.
            print(f'\n# quadratic intrastate parameters', file=fop)
            for i in range(nvib):
                k = Hq[i,i]
                print(f'  KD{state:02g}_{state:02g}_{i+1:03g}_{i+1:03g}  =     {k:12.8f}  ,  ev', file=fop)
            print(f'\n# bilinear intrastate parameters', file=fop)
            for i in range(nvib):
                for j in range(i):
                    k = Hq[i,j]
                    print(f'  KD{state:02g}_{state:02g}_{j+1:03g}_{i+1:03g}  =     {k:12.8f}  ,  ev', file=fop)


    # Done with parameters

    # Print hamiltonial elements
    section='''
end-parameter-section

HAMILTONIAN-SECTION
modes | el'''
    print(section, file=fop)

    # Define DoFs
    for i in range(nvib):
        print(f'modes | v{i+1:03g}', file=fop)
    print('', file=fop)

    # Kinetic
    for i in range(nvib):
        print(f'1.0*w{i+1:03g}    |{i+2:<3}   KE', file=fop)
    print('', file=fop)

    # Potential
    # Second order term (q²)
    for i in range(nvib):
        print(f'0.5*w{i+1:03g}    |{i+2:<3}   q^2', file=fop)
    print('', file=fop)
    # Consant term
    for state,null in enumerate(Ener,1):
        print(f'delta{state:02g}_{state:02g}        |1 S{state}&{state}', file=fop)
    print('', file=fop)
    # First order term (q)
    for state,gradq in enumerate(Gradq,1):
        for i,g in enumerate(gradq):
            if np.abs(g) > zero:
                print(f'K{state:02g}_{state:02g}_{i+1:03g}    |1 S{state}&{state}  |{i+2:<3}  q', file=fop)
        print('', file=fop)

    for state,Hq in enumerate(Hqs,1):
        # If VH, include second order contributions
        if Hq is not None:
            print('', file=fop)
            for i in range(nvib):
                print(f'KD{state:02g}_{state:02g}_{i+1:03g}_{i+1:03g}    |1 S{state}&{state}  |{i+2:<3}  q^2', file=fop)
            print('', file=fop)
            for i in range(nvib):
                for j in range(i):
                    print(f'2.0*KD{state:02g}_{state:02g}_{j+1:03g}_{i+1:03g}    |1 S{state}&{state}  |{j+2:<3}  q |{i+2:<3} q', file=fop)

    # Second order terms (ES freqs)
    # End file
    section='''
end-hamiltonian-section

end-operator
    '''
    print(section, file=fop)

def write_input(finp,nvib,integrator='mctdh'):
    section='''#######################################################################
###           Input Quantics
#######################################################################

RUN-SECTION

 name = qd
 propagate
 auto = twice
 steps
 psi = double
 gridpop

 tfinal= 300.0
 tout=   0.5
 tpsi=   5.0

end-run-section

OPERATOR-SECTION
 opname = oper
end-operator-section

sbasis-section                                                                  '''
    print(section, file=finp)

    for i in range(nvib):
        print(f'    v{i+1:03g}   =  10', file=finp)
    print('end-sbasis-section', file=finp)
    print('', file=finp)

    print('pbasis-section', file=finp)
    for i in range(nvib):
        print(f'    v{i+1:03g}   HO   20   0.0   1.0   1.0', file=finp)
    # Set number of states to 2 (even if only one is used to avoid issues)
    print('    el     el   2', file=finp)
    print('end-pbasis-section', file=finp)
    print('', file=finp)

    if integrator == 'mctdh':
        section='''INTEGRATOR-SECTION
 CMF/var =   0.05 ,  1.0E-5
 BS/spf  =     9 ,   1.0E-5
 SIL/A   =    50 ,   1.0E-5
end-integrator-section
'''
    else:
        section='''INTEGRATOR-SECTION
 VMF
 RK8 = 1.0d-7, 1.0d-4
end-integrator-section

'''

    section+='''INIT_WF-SECTION
 build
 init_state =    1                                                               '''
    print(section, file=finp)

    for i in range(nvib):
        print(f'    v{i+1:03g}   HO  0.0   0.0   1.0', file=finp)
    section='''end-build
end-init_wf-section

end-input
    '''
    print(section, file=finp)


if __name__ == '__main__':

    # Manage command line args
    args = sys.argv
    if len(args) == 1:
        print(f'Usage: \n {sys.argv[0]} fccout_fname(s) [-o oper_fname]')
        sys.exit()

    # remove app name from args
    args.pop(0)

    # Get fname if given
    fname_out = 'oper.op'
    if '-o' in args:
        i = args.index('-o')
        args.pop(i)
        fname_out = args.pop(i)

    # The remaining files shuld be the outputs
    fnames = args
    for fname in fnames:
        try:
            f = open(fname)
            f.close()
        except:
            print('Output from FCclasses not found: {fname}\n')
            print(f'Usage: \n {sys.argv[0]} fccout_fname(s) [-o oper_fname]')
            sys.exit()

    # Get data
    omega, Ener, Gradq, Hqs = [], [], [], []
    for fname in fnames:
        w,e,g,h = get_state_params(fname)
        omega = w
        Ener.append(e)
        Gradq.append(g)
        Hqs.append(h)

    # Generate operator
    fop = open(fname_out,'w')
    write_op(fop, omega, Ener, Gradq, Hqs)
    fop.close()

    # Generate tentative input
    nvib = len(omega)
    finp = open('qd.inp','w')
    write_input(finp,nvib,integrator='ml')
    finp.close()

