import numpy as np
import math

#Ex = input("Ex(Gpa):")

def D_to_R(teta):
    global R_teta
    R_teta = ((math.pi)*teta)/180

teta1 = 0
teta2 = 30
teta3 = -45

Q1 = np.matrix('181.8 2.897 0; 2.897 10.35 0; 0 0 7.17')

def Q_(Q, teta):
    tetat = D_to_R(teta)
    
    Q_11 = Q[0,0]*math.cos(R_teta)**4 + Q[1,1]*math.sin(R_teta)**4 + 2*(Q[0,1] + 2*Q[2,2])*(math.sin(R_teta)**2)*(math.cos(R_teta)**2)
    Q_12 = ( Q[0,0] + Q[1,1] -4*Q[2,2])*(math.sin(R_teta)**2)*(math.cos(R_teta)**2) + Q[0,1]*(math.cos(R_teta)**4 + math.sin(R_teta)**4)
    Q_16 = ( Q[0,0] - Q[0,1] -2*Q[2,2] )*(math.cos(R_teta)**3)*(math.sin(R_teta)) - ( Q[1,1] - Q[0,1] - 2*Q[2,2] )*(math.sin(R_teta)**3)*(math.cos(R_teta))

    Q_21 = Q_12
    Q_22 = (Q[0,0]*math.sin(R_teta)**4) + (Q[1,1]*math.cos(R_teta)**4) + 2*( Q[0,1] + 2*Q[2,2] )*(math.sin(R_teta)**2)*(math.cos(R_teta)**2)
    Q_26 = ( Q[0,0] - Q[0,1] -2*Q[2,2] )*((math.cos(R_teta))*(math.sin(R_teta)**3)) - ( Q[1,1] - Q[0,1] -2*Q[2,2] )*(math.cos(R_teta)**3)*(math.sin(R_teta))
    Q_66 = (Q[0,0] + Q[1,1] - 2*Q[0,1] -2*Q[2,2])*(math.sin(R_teta)**2)*(math.cos(R_teta)**2) + Q[2,2]*(math.sin(R_teta)**4 + math.cos(R_teta)**4)
    globals().update(locals())

Q_(Q1, teta3)

print(Q_16)

