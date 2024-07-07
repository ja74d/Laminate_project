import numpy as np
import math

#Ex = input("Ex(Gpa):")

def D_to_R(teta):
    global R_teta
    R_teta = ((math.pi)*teta)/180
    return R_teta

# Stacking Sequence
SS = [0, 30,-45]

Q1 = np.matrix('181.8 2.897 0; 2.897 10.35 0; 0 0 7.17')*1e+09

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
    TR[1,2] = -2*(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,0] = -(math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,1] = (math.cos(D_to_R(deg)))*(math.sin(D_to_R(deg)))
    TR[2,2] = math.cos(D_to_R(deg))**2 - math.sin(D_to_R(deg))**2
    globals().update(locals())
    return TR


#Stiffness Matrix of the laminate

#k = np.block([ [A, B], [B, A] ])

k = np.zeros((6, 6))

k[0,0] = A[0,0]
k[0, 1] = k[1, 0] = A[0, 1]
k[0, 2] = k[2, 0] = A[0, 2]
k[0, 3] = k[3, 0] = B[0, 0]
k[0, 4] = k[4, 0] = B[0, 1]
k[0, 5] = k[5, 0] = B[0, 2]

k[1, 1] = A[1, 1]
k[1, 2] = k[2, 1] = A[1, 2]
k[1, 3] = k[3, 1] = B[1, 0]
k[1, 4] = k[4, 1] = B[1, 1]
k[1, 5] = k[5, 1] = B[1, 2]

k[2, 2] = A[2, 2]
k[2, 3] = k[3, 2] = B[2, 0]
k[2, 4] = k[4, 2] = B[2, 1]
k[2, 5] = k[5, 2] = B[2, 2]

k[3, 3] = D[0, 0]
k[3, 4] = k[4, 3] = D[0, 1]
k[3, 5] = k[5, 3] = D[0, 2] 

k[4, 4] = D[1, 1]
k[4, 5] = k[5, 4] = D[1, 2]

k[5, 5] = D[2, 2]

k_prime = np.linalg.inv(k)

N = np.matrix('1000; 1000; 0; 0; 0; 0')

#mid-plane strain and curve

e0k = k_prime @ N

e0 = e0k[0:3, 0:1]
kapa = e0k[3:6, 0:1]

repeated_list = [num for num in h[1:-1] for _ in range(2)]
h = [h[0]] + repeated_list + [h[-1]]

SS = [num for num in SS for _ in range(2)]
print(SS)

a = 0
strains = []
stresses = []
for z in h:
    strain = e0 + (z)*kapa
    stress = Q_(Q1, SS[a]) @ strain
    a += 1
    strains.append(strain)
    stresses.append(stress)

#locals
local_stresses = []
x = 0
for r in stresses:
    local_stress = TR_M(SS[x]) @ r
    x += 1
    local_stresses.append(local_stress)
print(local_stresses)