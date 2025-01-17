#!/usr/bin/env python
import numpy as np
import sys

# Constants and conversions
evtown = 8065.5446811132
autoev = 27.2113961
autown = 2.1947463068e5

if len(sys.argv) == 2:
    fname = sys.argv[1]
    fname_out = 'oper.op'
elif len(sys.argv) == 3:
    fname = sys.argv[1]
    fname_out = sys.argv[2]
else:
    print(f'Usage: \n {sys.argv[0]} fccout_fname [oper_fname]')
    sys.exit()

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
    # Find location of Duschinsky file
    dusch_file = None
    for line in f:
        if 'Printing Duschinski matrix to' in line:
            dusch_file = line.strip().split()[-1].replace("'","")
else:
    print('ERROR: only VG and VH models supported')
    sys.exit()
f.close()

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
omega *= autoev
for i,w in enumerate(omega):
    print(f'  w{i+1:03g}  =     {w:12.8f}  ,  ev', file=fop)

# Add energy of the state (=0)
section='''
# vertical energies and off diag coupling const (if present)
    delta01_01 =        0.00000000  ,   ev

# linear intrastate parameters'''
print(section, file=fop)


# Print grad in eV
gradq *= autoev
for i,g in enumerate(gradq):
    print(f'  KD01_01_{i+1:03g}  =     {g:12.8f}  ,  ev', file=fop)

# If VH, the modify ES forces
if model == 'VH':
    Hq *= autoev / 2.
    print(f'\n# quadratic intrastate parameters', file=fop)
    for i in range(nvib):
        k = Hq[i,i] 
        print(f'  KD01_01_{i+1:03g}_{i+1:03g}  =     {k:12.8f}  ,  ev', file=fop)
    print(f'\n# bilinear intrastate parameters', file=fop)
    for i in range(nvib):
        for j in range(i):
            k = Hq[i,j] 
            print(f'  KD01_01_{j+1:03g}_{i+1:03g}  =     {k:12.8f}  ,  ev', file=fop)


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
print('delta01_01        |1 S1&1', file=fop)
print('', file=fop)
# First order term (q)
for i in range(nvib):
    print(f'KD01_01_{i+1:03g}    |1 S1&1  |{i+2:<3}  q', file=fop)
# If VH, include second order contributions 
if model == 'VH':
    print('', file=fop)
    for i in range(nvib):
        print(f'KD01_01_{i+1:03g}_{i+1:03g}    |1 S1&1  |{i+2:<3}  q^2', file=fop)
    print('', file=fop)
    for i in range(nvib):
        for j in range(i):
            print(f'2.0*KD01_01_{j+1:03g}_{i+1:03g}    |1 S1&1  |{j+2:<3}  q |{i+2:<3} q', file=fop)

# Second order terms (ES freqs)
# End file
section='''
end-hamiltonian-section

end-operator
'''
print(section, file=fop)

fop.close()

