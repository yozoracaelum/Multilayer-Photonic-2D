#Multilayer 1D
import matplotlib.pyplot as plt
import numpy as np

def E1(z):
    A1 = 0
    B1 = 0
    w = 0
    n1 = 0
    k1 = w*n1
    result = A1*np.exp(1j*k1*z) + B1*np.exp(-1j*k1*z)
    return result

def E2(z):
    C1 = 0
    D1 = 0
    w = 0
    n2 = 0
    k2 = w*n2
    d1 = 0
    result = C1*np.exp(1j*k2*(z-d1)) + D1*np.exp(-1j*k2*(z-d1))
    return result

n1 = float(input("n1: "))
n2 = float(input("n2: "))
lamda = np.arange(1,1500,0.001)
w = (2*np.pi)/lamda
#k0 = 0
k1 = w*n1
k2 = w*n2
d1 = 10
d2 = 10

M11 = np.exp(1j*k1*d1)*(np.cos(k2*d2)+0.5j*((k2/k1)+(k1/k2))*np.sin(k2*d2))

'''
M12 = np.exp(-1j*k1*d1)*(0.5j*((k2/k1)-(k1/k2))*np.sin(k2*d2))
M21 = np.exp(1j*k1*d1)*(-0.5j*((k2/k1)-(k1/k2))*np.sin(k2*d2))
M22 = np.exp(-1j*k1*d1)*(np.cos(k2*d2)-0.5j*((k2/k1)+(k1/k2))*np.sin(k2*d2))

t = 2*1j*k1*np.exp(-1j*k2*d1)*(1/(-M21+k1*k2*M12+1j*(k2*M11+k1*M22)))
'''
'''
M11 = np.exp(1j*k1*d1)*(np.cos(k2*d2)+0.5j*((k2/k1)+(k1/k2))*np.sin(k2*d2))
M12 = np.exp(-1j*k1*d1)*(0.5j*((k2/k1)-(k1/k2))*np.sin(k2*d2))
M21 = np.exp(1j*k1*d1)*(-0.5j*((k2/k1)-(k1/k2))*np.sin(k2*d2))
M22 = np.exp(-1j*k1*d1)*(np.cos(k2*d2)-0.5j*((k2/k1)+(k1/k2))*np.sin(k2*d2))

M = [
    [M11,M12],
    [M21,M22]
    ]

term1 = [
        [0.5*(1+(k0/k1)),0.5*(1-(k0/k1))],
        [0.5*(1-(k0/k1)),0.5*(1+(k0/k1))]
        ]
term2 = [
        [0.5*(1+(k2/k1)),0.5*(1-(k2/k1))],
        [0.5*(1-(k2/k1)),0.5*(1+(k2/k1))]
        ]
term3 = [
        [0.5*(1+(k2/k0)),0.5*(1-(k2/k0))],
        [0.5*(1-(k2/k0)),0.5*(1+(k2/k0))]
        ]

Mt1 = np.dot(term1,(np.dot(M,np.linalg.inv(term2))))
Mt = np.dot(Mt1,term3)
'''

T = (1/M11)**2
print(T)
print(len(lamda))
print(len(T))
#print(Mt)
plt.plot(lamda,T)
plt.show()
