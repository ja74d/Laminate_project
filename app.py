import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

E1 = 181e+09
E2 = 10.3e+09
v12 = 0.25
v21 = ((v12)*E2)/(E1)
G12 = 7.17e+09

vm = 0.3
vf = 0.3

Em = 3.4e+09
Ef = 230e+09
Gf = Ef/(2*(1+vf))
Gm = Em/(2*(1+vm))


Vf = 0.7
Vm = 0.3
def HT():
    #Halpine Tsi
    ksi = 1 + (40*(Vf**10))
    eta = (Ef-Em)/(Ef + (ksi*Em))
    E2 = ((Em)*(1+ksi*eta*Vf))/(1-eta*Vf)
    print(E2)

    #G12
    eta = (Gf-Gm)/(Gf + (ksi*Gm))
    G12 = ((Gm)*(1+ksi*eta*Vf))/(1-eta*Vf)
    #print(G12)

Q = np.zeros((3, 3))
S = np.zeros((3, 3))

S[0, 0] = 1/E1
S[0, 1] = -v12/E1
S[0, 2] = 0

S[1, 0] = S[0, 1]
S[1, 1] = 1/E2
S[1, 2] = 0

S[2, 0] = 0
S[2, 1] = 0
S[2, 2] = 1/G12

Q1 = np.linalg.inv(S)


def D_to_R(teta):
    global R_teta
    R_teta = ((math.pi)*teta)/180
    return R_teta

# Stacking Sequence
SS = [0, 30, -45]
SST = [num for num in SS for _ in range(2)]

#Q1 = np.matrix('181.8 2.897 0; 2.897 10.35 0; 0 0 7.17')*1e+09

#Reduced Stiffnes Matrix
def Q_(Q, teta):
    tetat = D_to_R(teta)
    
    Q_11 = Q[0,0]*math.cos(R_teta)**4 + Q[1,1]*math.sin(R_teta)**4 + 2*(Q[0,1] + 2*Q[2,2])*(math.sin(R_teta)**2)*(math.cos(R_teta)**2)
    Q_12 = ( Q[0,0] + Q[1,1] -4*Q[2,2])*(math.sin(R_teta)**2)*(math.cos(R_teta)**2) + Q[0,1]*(math.cos(R_teta)**4 + math.sin(R_teta)**4)
    Q_16 = ( Q[0,0] - Q[0,1] -2*Q[2,2] )*(math.cos(R_teta)**3)*(math.sin(R_teta)) - ( Q[1,1] - Q[0,1] - 2*Q[2,2] )*(math.sin(R_teta)**3)*(math.cos(R_teta))

    Q_21 = Q_12
    Q_22 = (Q[0,0]*math.sin(R_teta)**4) + (Q[1,1]*math.cos(R_teta)**4) + 2*( Q[0,1] + 2*Q[2,2] )*(math.sin(R_teta)**2)*(math.cos(R_teta)**2)
    Q_26 = ( Q[0,0] - Q[0,1] -2*Q[2,2] )*((math.cos(R_teta))*(math.sin(R_teta)**3)) - ( Q[1,1] - Q[0,1] -2*Q[2,2] )*(math.cos(R_teta)**3)*(math.sin(R_teta))

    Q_61 = Q_16
    Q_62 = Q_26
    Q_66 = (Q[0,0] + Q[1,1] - 2*Q[0,1] -2*Q[2,2])*(math.sin(R_teta)**2)*(math.cos(R_teta)**2) + Q[2,2]*(math.sin(R_teta)**4 + math.cos(R_teta)**4)
    Q__ = np.matrix([[Q_11 ,Q_12 ,Q_16], [Q_21, Q_22, Q_26], [Q_16, Q_26, Q_66]])
    return Q__
    globals().update(locals())
    

#Q_(Q1, SS[1])

#locations of the ply surfaces

#h = [-0.0075, -0.0025, 0.0025, 0.0075]

# number of layers
nl = len(SS)

#ply thickness
h_ = 0.005

midplane = nl*h_
h0 = ((-1)*midplane/2)
h = [h0]

for n in range(0, nl):
    hi = h[n] + h_
    h.append(hi)

#for er in range(len(SS)):
#    print(f'teta = {SS[er]}:\n', Q_(Q1, SS[er]))

# A Matrix

A = np.zeros((3,3))

for i in range(0,len(SS)):
    A += Q_(Q1, SS[i])*(h[i+1] - h[i])

A_prime = np.linalg.inv(A)



#print('A:', A)
#print()
# B Matrix

B = np.zeros((3, 3))

for j in range(0, len(SS)):
    B += 0.5*(Q_(Q1, SS[j])*((h[j+1])**2 - (h[j])**2))

B_prime = np.linalg.inv(B)

#print('B:', B)
#print()
# D Matrix

D = np.zeros((3, 3))

for ka in range(0, len(SS)):
    D += (1/3)*(Q_(Q1, SS[ka])*((h[ka+1])**3 - (h[ka])**3))

D_prime = np.linalg.inv(D)

#print('D:', D)
#print()

#in-plane Engineering constants

Ex = 1/(A_prime[0,0]*midplane)

Ey = 1/(A_prime[1,1]*midplane)

Gxy = 1/(A_prime[2,2]*midplane)

Vxy = -( (A_prime[0,1])/(A_prime[0,0]) )

Vyx = -( (A_prime[0,1])/(A_prime[1,1]) )

#Flexural Engineering constants

Exf = 12/( (midplane**3)*(D_prime[0,0]) )

Eyf = 12/((midplane**3)*(D_prime[1,1]))

Gxyf = 12/((midplane**3)*(D_prime[2,2]))

Vxyf = -( (D_prime[0,1])/(D_prime[0,0]) )

Vyxf = -( (D_prime[0,1])/(D_prime[1,1]) )




#the matrix that moves global stress to local stress

TR = np.zeros((3, 3))

def TR_M(deg):
    TR[0,0] = (math.cos(D_to_R(deg)))**2
    TR[0,1] = (math.sin(D_to_R(deg)))**2
    TR[0,2] = 2*(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[1,0] = TR[0,1]
    TR[1,1] = (math.cos(D_to_R(deg)))**2 
    TR[1,2] = -2*(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,0] = -(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,1] = (math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,2] = math.cos(D_to_R(deg))**2 - math.sin(D_to_R(deg))**2
    globals().update(locals())
    return TR


#Stiffness Matrix of the laminate

def ABD(a11, b11, c11):
    global k
    global k_prime
    # Compute the 2x2 block matrix
    k = np.block([[a11, b11], [b11, c11]])

    # Flatten the matrix into a 1D array
    ka = k.flatten()

    # Convert the flattened array back to a 6x6 matrix
    k = ka.reshape(6, 6)

    # Compute the inverse of k
    k_prime = np.linalg.inv(k)

    return k, k_prime


ABD(A, B, D)
#print(k)

N = np.matrix('1000; 1000; 0; 0; 0; 0')

#mid-plane strain and curve


 # Extend the h list
repeated_list = [num for num in h[1:-1] for _ in range(2)]
h_extended = [h[0]] + repeated_list + [h[-1]]

# Extend the SS list
SST = [num for num in SS for _ in range(2)]


def ek(k):
    global local_strains, local_stresses, e0, kapa, strains, stresses
    # Calculate e0 and kapa
    e0k = np.linalg.inv(k) @ N
    e0 = e0k[0:3, 0:1]
    kapa = e0k[3:6, 0:1]

   
    # Calculate strains and stresses
    strains = []
    stresses = []
    for a, z in enumerate(h_extended):
        strain = e0 + z * kapa
        stress = Q_(Q1, SST[a]) @ strain
        strains.append(strain)
        stresses.append(stress)

    # Calculate local stresses
    local_stresses = [TR_M(SST[x]) @ stresses[x] for x in range(len(SST))]

    # Calculate local strains
    local_strains = []
    for y, f in enumerate(SST):
        strn = np.matrix('1, 0, 0; 0, 1, 0; 0, 0, 0.5') @ strains[y]
        local_strain = TR_M(SST[y]) @ strn
        local_strains.append(local_strain)

    return local_strains, local_stresses
ek(k)
#print('e0:\n', e0)
#print()
#print('kapa:\n', kapa)

#tables for stress and strain in top and bottom for each ply
name = [f'{SS[0]} top', f'{SS[0]} bottom', f'{SS[1]} top', f'{SS[1]} bottom', f'{SS[2]} top', f'{SS[2]} bottom']

data = {name: mat.flatten().tolist()[0] for name, mat in zip(name, local_strains)}
df = pd.DataFrame(data, index=['epsilon 1', 'epsilon 2', 'gama'])

# Print the DataFrame
#print(df)

#T and C (hygrothermal)

delta_T = -75

alpha_local = np.matrix('0.2e-07; 0.225e-04; 0')

alpha_global = []
strain_T = []

for e in range(len(SS)):
    alpha = TR_M(SS[e]) @ alpha_local
    alpha_global.append(alpha)
    stn = alpha_global[e]*delta_T
    strain_T.append(stn)
#print(strain_T)

repeated_list3 = [num for num in strain_T[1:-1] for _ in range(2)]
strain_T = [strain_T[0]] + repeated_list + [strain_T[-1]]


Ntc = np.zeros((3, 1))
Mtc = np.zeros((3, 1))

#ply thickness
h_ = 0.005

midplane = nl*h_
h0 = ((-1)*midplane/2)
h = [h0]

for n in range(0, nl):
    hi = h[n] + h_
    h.append(hi)


for u in range(0, len(SS)):
    Ntc += delta_T*(Q_(Q1, SS[u])@alpha_global[u])*(h[u+1] - h[u])

for v in range(0, len(SS)):
    Mtc += (0.5*delta_T)*(Q_(Q1, SS[v])@alpha_global[v])*((h[v+1])**2 - (h[v])**2)

MN = np.zeros((6, 1))
MN[0, 0] = Ntc[0, 0]
MN[1, 0] = Ntc[1, 0]
MN[2, 0] = Ntc[2, 0]

MN[3, 0] = Mtc[0, 0]
MN[4, 0] = Mtc[1, 0]
MN[5, 0] = Mtc[2, 0]

repeated_list = [num for num in h[1:-1] for _ in range(2)]
h = [h[0]] + repeated_list + [h[-1]]

e0ktc = k_prime@MN

e0_tc = e0ktc[0:3, 0:1]
kapa_tc = e0ktc[3:6, 0:1]


strains_tc = []

for z in h:
    strain = e0_tc + (z)*kapa_tc
    strains_tc.append(strain)

#print(strains_tc)

#print(strain_T)

# Mechanical Strains
M_strains = []
for f in range(0, len(SST)):
    M_strain = strains_tc[f] - strain_T[f]
    M_strains.append(M_strain)
#print(M_strains)

# M stresses
repeated_list = [num for num in h[1:-1] for _ in range(2)]
h = [h[0]] + repeated_list + [h[-1]]

SST = [num for num in SS for _ in range(2)]

stresses = []

a = 0
for z in SST:
    stress = Q_(Q1, SST[a]) @ M_strains[a]
    stresses.append(stress)
    a += 1
#print(stresses)


# Tsi-Wu Failure theory

F1t = 1500e+06
F1c = 1500e+06
F2t = 40e+06
F2c = 246e+06
F6 = 68e+06

f1 = (1/(F1t)) - (1/(F1c))
f11 = 1/((F1t)*(F1c))

f2 = (1/(F2t)) - (1/(F2c))
f22 = 1/((F2t)*(F2c))

f66 = 1/((F6)*(F6))

f12 = (-0.5) * (1/(F1c*F1t*F2t*F2c))**0.5
#f12 = -1/(2*((F1t)*(F1t)))
#f12 = -3.36032e-18

#maximum stresses
max_stress = []
for strs in stresses:
    sigmax = float(strs[0])
    sigmay = float(strs[1])
    tau = float(strs[2])
    simgamxn1 = ((sigmax+sigmay)/2) + (((sigmax-sigmay)/2)**2+(tau)**2)**0.5
    sigmamxn2 = ((sigmax+sigmay)/2) - (((sigmax-sigmay)/2)**2+(tau)**2)**0.5
    taumax = (((sigmax-sigmay)/2)**2+(tau)**2)**0.5
    sigmamxn = [simgamxn1, sigmamxn2, taumax]
    max_stress.append(sigmamxn)
#print(max_stress)

#maximum strains
max_strain = []
for strn in strains:
    epsx = float(strn[0])
    epsy = float(strn[1])
    gama = float(strn[2])
    epsmxn1 = ((epsx+epsy)/2) + (((epsx-epsy)/2)**2+(gama)**2)**0.5
    epsmxn2 = ((epsx+epsy)/2) - (((epsx-epsy)/2)**2+(gama)**2)**0.5
    gamamax = (((epsx-epsy)/2)**2+(gama)**2)**0.5
    epsmxn = [epsmxn1, epsmxn2, gamamax]
    max_strain.append(epsmxn)
#print(max_strain)


def Tsiwu(local_stresses):
    global Tsi_wu, Tsi_wu_stress
    Tsi_wu = []
    Tsi_wu_stress = []
    
    for st in local_stresses:
        sigma1 = float(st[0])
        sigma2 = float(st[1])
        tau = float(st[2])

        # Calculate TW
        TW = ((f1)*(sigma1)) + ((f2)*(sigma2)) + ((f11)*(sigma1)**2) + ((f22)*(sigma2)**2) + ((f66)*(tau)**2) + ((2)*(f12)*(sigma1)*(sigma2))

        # Calculate SR
        a = (f11*((sigma1)**2)) + (f22*((sigma2)**2)) + (f66*((tau)**2)) + (2*f12*sigma1*sigma2)
        b = (f1*sigma1) + (f2*sigma2)
        c = -1 

        SR = np.roots([a, b, c])
        
        for root in SR:
            if root > 0:
                Tsi_wu.append(root)
                Tsi_wu_stress.append(root / (len(SS) * h_))

    return Tsi_wu, Tsi_wu_stress

Tsiwu(local_stresses)

#Tsi-hill
def Tsihill(local_stresses):
    global Tsi_hill, Tsi_hill_stress
    Tsi_hill = []
    Tsi_hill_stress = []
    aa = 0
    for ts in local_stresses:
        sigma1 = local_stresses[aa][0]
        sigma1 = float(sigma1)
        sigma2 = local_stresses[aa][1]
        sigma2 = float(sigma2)
        tau = local_stresses[aa][2]
        tau = float(tau)
        aa += 1

        TH = (1/((sigma1/F1t)**2+(sigma2/F2t)**2+(tau/F6)**2-((sigma1*sigma2)/(F1t**2))))**0.5
        Tsi_hill.append(TH)
    for g in Tsi_hill:
        vb = g/((len(SS))*(h_))
        Tsi_hill_stress.append(vb)

    return Tsi_hill, Tsi_hill_stress
    #print(Tsi_hill)
    #print(Tsi_hill_stress)
    #print(Tsi_hill_stress)

Tsihill(local_stresses)

first = Tsi_hill_stress.index(min(Tsi_hill_stress))

first_layar_failure = SS[int(first/2)]


first_ply_stress = Tsi_wu_stress[first]

first_ply_strain = min(Tsi_wu) * e0[0]

# now how the fuck i`m spouse to set the stiffness of this layer to zero?

new_stiffness = []
for lay in SS:
    if lay == first_layar_failure:
        #print(np.zeros((3, 3)))
        new_stiffness.append(np.zeros((3, 3)))
    else:
        #print(Q_(Q1, lay))
        new_stiffness.append(Q_(Q1, lay))
#print(new_stiffness)

midplane = nl*h_
h0 = ((-1)*midplane/2)
h = [h0]

for n in range(0, nl):
    hi = h[n] + h_
    h.append(hi)

#New A B D

def new_A(SS, stiffness, h):
    global NA
    global NA_prime
    # Initialize NA as a 3x3 zero matrix
    NA = np.zeros((3, 3))

    # Calculate NA
    for ik in range(0, len(SS)):
        NA += stiffness[ik] * (h[ik+1] - h[ik])

    # Calculate the inverse of NA
    NA_prime = np.linalg.inv(NA)

    return NA, NA_prime

new_A(SS, new_stiffness, h)
#print(NA)

def new_B(SS, stiffness, h):
    global NB
    global NB_prime
    # Initialize NB as a 3x3 zero matrix
    NB = np.zeros((3, 3))

    # Calculate NB
    for nj in range(0, len(SS)):
        NB += 0.5 * stiffness[nj] * ((h[nj+1]) ** 2 - (h[nj]) ** 2)

    # Calculate the inverse of NB
    NB_prime = np.linalg.inv(NB)

    return NB, NB_prime
new_B(SS, new_stiffness, h)
#print(NB)

def new_D(SS, stiffness, h):
    global ND
    global ND_prime
    # Initialize ND as a 3x3 zero matrix
    ND = np.zeros((3, 3))

    # Calculate ND
    for nk in range(0, len(SS)):
        ND += (1/3) * stiffness[nk] * ((h[nk+1]) ** 3 - (h[nk]) ** 3)

    # Calculate the inverse of ND
    ND_prime = np.linalg.inv(ND)

    return ND, ND_prime
new_D(SS, new_stiffness, h)

ABD(NA, NB, ND)

#print(k)

ek(k)


Tsiwu(local_stresses)

#print(first)

del Tsi_wu_stress[first]
if first == (len(SST) - 1):
    del Tsi_wu_stress[first-1]
else:
    del Tsi_wu_stress[first]

del Tsi_wu[first]
if first == (len(SST) - 1):
    del Tsi_wu[first-1]
else:
    del Tsi_wu[first]

#print(local_stresses)

second = Tsi_wu_stress.index(min(Tsi_wu_stress))

second_layar_failure = SS[int(second/2)]

second_ply_stress = Tsi_wu_stress[second]

second_ply_strain = min(Tsi_wu) * e0[0]



xd = [0, float(first_ply_strain), float(second_ply_strain)]
yd = [0, float(first_ply_stress), float(second_ply_stress)]

plt.plot(xd, yd, marker='o')

#plt.show()

del SS[int(second/2)]
del SS[int(second/2)]

last_ply = SS[0]
#print(last_ply)