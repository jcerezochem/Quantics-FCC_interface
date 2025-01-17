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
else:
    print(f'Usage: \n {sys.argv[0]} oper_fname')
    sys.exit()

f = open(fname)

#############
# PARAMETERS
#############
# Find parameter section
for line in f:
    line=line.split('#')[0]
    line=line.strip()
    if 'PARAMETER-SECTION' in line.upper():
        break

# Read parameters
parameters = {}
for line in f:
    line=line.split('#')[0]
    line=line.strip()
    if 'END' in line.upper():
        break
    if not line:
        continue
    label,data = line.split('=')
    label=label.replace(' ','')
    if ',' in data:
        value,units = data.split(',')
    else:
        value,units = [data, None]
    value = float(value)
    units = units.replace(' ','')
    parameters[label] = [value,units]
    
##############
# HAMILTONIAN
##############
for line in f:
    line=line.split('#')[0]
    line=line.strip()
    if 'HAMILTONIAN-SECTION' in line.upper():
        break

# Read Hamiltonian
# Frist get modes and matrix elements ids
modes=[]
matelems=[]
i_el = 0
for line in f:
    line=line.split('#')[0]
    line=line.strip()
    if 'END' in line.upper():
        break
    if not line:
        continue
    if 'MODES' in line.upper():
        modes += line.split('|')[1:]
    elif '|' in line:
        data = line.split('|')
        data.pop(0)
        for entry in data:
            i, label = entry.split()
            i = int(i)
            if 'S' in label:
                i_el = i
                if label not in matelems:
                    matelems.append(label)
                break

# Rewind file
f.close()
f = open(fname)

# Manage number of modes
if i_el > 0:
    nvib = len(modes) - 1
else:
    nvib = len(modes)

# Manage diagonal and off-diagonal el matrix
omega = np.zeros(nvib)
Hq = {}
gradq = {}
v = {}
for matelem in matelems:
    Hq[matelem] = np.zeros([nvib,nvib])
    gradq[matelem] = np.zeros(nvib)
    v[matelem] = 0

# Map modes to matrix indices
map_modes = {}
k = 0
for i,j in enumerate(modes,1):
    if i == i_el:
        continue
    map_modes[i] = k
    k+=1
        
# Read Hamiltonian
for line in f:
    line=line.split('#')[0]
    line=line.strip()
    if 'HAMILTONIAN-SECTION' in line.upper():
        break

for line in f:
    line=line.split('#')[0]
    line=line.strip()
    if 'END' in line.upper():
        break
    if not line:
        continue
    if 'MODES' in line.upper():
        pass
    elif '|' in line:
        data = line.split('|')
        value = data[0].replace(' ','')
        mult = 1.
        if '*' in value:
            mult,value = value.split('*')
            mult=float(mult)
        elif value[0] == '-':
            value=value[1:]
            mult=-1.
        value=mult*parameters[value][0]
        # Get the matrix element
        data.pop(0)
        matelem = 'ALL'
        hamilton_term = []
        coords = []
        for entry in data:
            i, label = entry.split()
            i = int(i)
            if i == i_el:
                matelem = label
            else:
                coords.append(map_modes[i])
                hamilton_term.append(label)
        if matelem == 'ALL':
            if hamilton_term == ['q^2']:
                i = coords[0]
                omega[i] += value / autoev * 2.
            elif hamilton_term == ['KE']:
                print('KE entry ignored:', line)
            else:
                print('Unexpected entry:', line)
        else:
            if hamilton_term == []:
                v[matelem] += value / autoev
            elif hamilton_term == ['q']:
                i = coords[0]
                gradq[matelem][i] = value / autoev 
            elif hamilton_term == ['q^2']:
                i = coords[0]
                Hq[matelem][i,i] += value / autoev * 2.
            elif hamilton_term == ['q','q']:
                i,j = coords
                Hq[matelem][i,j] += value / autoev 
                Hq[matelem][j,i] += value / autoev 
            else:
                print('Unknown entry ignored:', line)
f.close()

# Change to Q
for matelem in matelems:
    i,j = [ int(m.replace('S','')) for m in matelem.split('&') ]
    if i == j:
        Hq[matelem] += np.diag(omega)
    gradq[matelem] *= np.sqrt(omega)
    for i in range(nvib):
        Hq[matelem][i,i] *= omega[i]
        for j in range(i):
            Hq[matelem][i,j] *= np.sqrt(omega[i]*omega[j])
            Hq[matelem][j,i] *= np.sqrt(omega[i]*omega[j])

# Write Freqs S0
f = open('Freq_S0.dat', 'w')
for w in omega:
    w *= autown
    print(w, file=f)
f.close()

# Compute ESs freqs and print
for matelem in matelems:
    i,j = [ int(m.replace('S','')) for m in matelem.split('&') ]
    if i != j:
        continue
    # Freqs
    omega_ES, J = np.linalg.eig(Hq[matelem])
    omega_ES = np.sqrt(omega_ES)
    f = open(f'Freq_S{i}.dat', 'w')
    for w in omega_ES:
        w *= autown
        print(w, file=f)
    f.close()
    # Grads 
    #gradq[matelem] = J.transpose() @ gradq[matelem] # keeep them in Q (for S0) coords!
    f = open(f'GradQ_S{i}.dat', 'w')
    for g in gradq[matelem]:
        print(g, file=f)
    f.close()
    # Dusch
    f = open(f'Dusch_S{i}.dat', 'w')
    for row in J: 
        for item in row:
            print(item, file=f)
    f.close()
    # Input
    input_text=f'''$$$
PROPERTY     =   OPA  ; OPA/EMI/ECD/CPL/RR/TPA/TPCD/MCD/IC
MODEL        =   VH   ; AS/ASF/AH/VG/VGF/VH
DIPOLE       =   FC   ; FC/HTi/HTf
DE           =   {v[matelem]}; In eV, Required with old state-files
TEMP         =   0.00 ; (temperature in K)
BROADFUN     =   GAU  ; GAU/LOR/VOI
HWHM         =   0.01 ; (broadening width in eV)
METHOD       =   TD   ; TI/TD
;VIBRATIONAL ANALYSIS
NVIB         =   {nvib}
NORMALMODES  =   IMPLICIT   ; COMPUTE/READ/IMPLICIT
;INPUT DATA FILES
ELDIP_FILE   =   eldip
FREQ1_FILE   =   Freq_S0.dat
FREQ2_FILE   =   Freq_S{i}.dat
DUSCH_FILE   =   Dusch_S{i}.dat
GRAD2_FILE   =   GradQ_S{i}.dat
GRAD2_COORD  =   NORMALMODE1
'''
    f = open(f'fcc_S{i}.inp', 'w')
    print(input_text, file=f)
    f.close()

f = open('eldip', 'w')
print('1. 0. 0.', file=f)
print('1. 0. 0.', file=f)
f.close()
 
