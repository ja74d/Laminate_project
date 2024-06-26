import numpy as np
import math

#Ex = input("Ex(Gpa):")

def D_to_R(teta):
    global R_teta
    R_teta = ((math.pi)*teta)/180
    return R_teta

# Stacking Sequence
SS = [0, 90, 0]

Q1 = np.matrix('181.8 2.897 0; 2.897 10.35 0; 0 0 7.17')

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


# A Matrix

A = np.zeros((3,3))

for i in range(0,3):
    A += Q_(Q1, SS[i])*(h[i+1] - h[i])

A_prime = np.linalg.inv(A)



#print(A)

# B Matrix

B = np.zeros((3, 3))

for j in range(0, 3):
    B += 0.5*(Q_(Q1, SS[j])*((h[j+1])**2 - (h[j])**2))

B_prime = np.linalg.inv(B)

#print(B)
# D Matrix

D = np.zeros((3, 3))

for k in range(0, 3):
    D += 0.3333*(Q_(Q1, SS[k])*((h[k+1])**3 - (h[k])**3))

D_prime = np.linalg.inv(D)

#print(D)

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

#T and C

Q_0 = Q_(Q1, 0)
#print(Q_0)

Q_90 = Q_(Q1, 90)
#print(Q_90)

alpha_local = np.matrix('0.2e-07; 0.225e-04; 0')

alpha_0 = alpha_local

alpha_90 = np.matrix('0.225e-04; 0.2e-07; 0')

delta_T = -75

#Failure of the laminate#


#the matrix that moves global stress to local stress

TR = np.zeros((3, 3))

def TR_M(deg):
    TR[0,0] = (math.cos(D_to_R(deg)))**2
    TR[0,1] = (math.sin(D_to_R(deg)))**2
    TR[0,2] = 2*(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[1,0] = TR[0,1]
    TR[1,1] = (math.cos(D_to_R(deg)))**2 
    TR[1,2] =-2*(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,0] = TR[0,2]
    TR[2,1] = TR[1,2]
    TR[2,2] = math.cos(D_to_R(deg))**2 - math.sin(D_to_R(deg))**2
    globals().update(locals())
    return TR


#local stress matrix for all layers in stacking sequence

# strains of each layer

N = np.matrix('1; 0; 0')

midplane_strains = A_prime*N

local_stresses = []
local_straines = []

for q in SS:
    TR_M(q)
    local_strain = TR*midplane_strains
    local_straines.append(local_strain)
    #print(local_strains)
    local_stress = TR*Q_(Q1, q)*midplane_strains
    local_stresses.append(local_stress)
    #print(local_Streses)



