import matplotlib.pyplot as plt
import numpy as np

n0 = 1
n1 = float(input("n1: "))
n2 = float(input("n2: "))
#l1 = float(input("l min: "))
#l2 = float(input("l max: "))
d = float(input("d: "))
N = int(input("N: "))
d1,d2 = d,d
m0,m1,m11,m00,m10,m10,m110,Mt = [],[],[],[],[],[],[],[]
lamda = np.arange(0.1,2,0.01)
w = [(2*np.pi)/i for i in lamda]
k0 = [i*n0 for i in w]
k1 = [i*n1 for i in w]
k2 = [i*n2 for i in w]
M11 = [np.exp(1j*k1[i]*d1)*(np.cos(k2[i]*d2)+0.5j*((k2[i]/k1[i])+(k1[i]/k2[i]))*np.sin(k2[i]*d2)) for i in range(len(lamda))]
M12 = [np.exp(-1j*k1[i]*d1)*(0.5j*((k2[i]/k1[i])-(k1[i]/k2[i]))*np.sin(k2[i]*d2)) for i in range(len(lamda))]
M21 = [np.exp(1j*k1[i]*d1)*(-0.5j*((k2[i]/k1[i])-(k1[i]/k2[i]))*np.sin(k2[i]*d2)) for i in range(len(lamda))]
M22 = [np.exp(-1j*k1[i]*d1)*(np.cos(k2[i]*d2)-0.5j*((k2[i]/k1[i])+(k1[i]/k2[i]))*np.sin(k2[i]*d2)) for i in range(len(lamda))]
for i in range(len(lamda)):
    m00.append(M11[i])
    m00.append(M12[i])
    m0.append(m00)
    m10.append(M21[i])
    m10.append(M22[i])
    m1.append(m10)
    m110.append(m0[i])
    m110.append(m1[i])
    m11.append(m110)
    m00,m10,m110 = [],[],[]
m11 = np.array(m11)
A = m11[0]
for i in range(len(lamda)):
    for j in range(N-1):
        A = A@m11[i]
    Mt.append(A)
Mt = np.array(Mt)
#print(Mt)
print(Mt[0][1][0])
r = [-(Mt[i][1][0]/Mt[i][1][1]) for i in range(len(lamda))]
t = [Mt[i][0][0]-((Mt[i][0][1]*Mt[i][1][0])/Mt[i][1][1]) for i in range(len(lamda))]
R = [np.abs(r[i])**2 for i in range(len(lamda))]
T = [np.abs(t[i])**2 for i in range(len(lamda))]
plt.title('n1: %.1f | n2: %.1f | d: %.2f micrometer | N: %d' %(n1,n2,d,N))
plt.plot(lamda,T,'b',label="Transmittance")
plt.plot(lamda,R,'r',label="Reflectance")
plt.xlabel("Lamda (micrometer)")
plt.ylabel("Intensity")
plt.grid()
plt.legend()
plt.show()

'''
term1 = [[
        [0.5*(1+(k0[i]/k1[i])),0.5*(1-(k0[i]/k1[i]))],
        [0.5*(1-(k0[i]/k1[i])),0.5*(1+(k0[i]/k1[i]))]
        ] for i in range(len(lamda))]
term2 = [[
        [0.5*(1+(k2[i]/k1[i])),0.5*(1-(k2[i]/k1[i]))],
        [0.5*(1-(k2[i]/k1[i])),0.5*(1+(k2[i]/k1[i]))]
        ] for i in range(len(lamda))]
term3 = [[
        [0.5*(1+(k2[i]/k0[i])),0.5*(1-(k2[i]/k0[i]))],
        [0.5*(1-(k2[i]/k0[i])),0.5*(1+(k2[i]/k0[i]))]
        ] for i in range(len(lamda))]
Minv = [np.linalg.inv(term2[i]) for i in range(len(lamda))]
Mt1 = [np.dot(term1[i],m11[i]) for i in range(len(lamda))]
Mt2 = [np.dot(Mt1[i],Minv[i]) for i in range(len(lamda))]
Mt3 = [np.dot(Mt2[i],term3[i]) for i in range(len(lamda))]
Mt = np.array(Mt3)
'''
