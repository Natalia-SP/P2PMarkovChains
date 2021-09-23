# Librerías
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def u():
  u = random.uniform(0,1000000)
  return u 

# Nodos infectados
Ni = 1

# Estado inicial

xn = 0
yn = 1
xi = Ni
yi = 1

# Parámetros

C=0.002
L = 2
M = 0.05
H = 0.005 
G = .3
N = 0.85
K = 1/10

Tsim = 0
pix = np.zeros([5000,6000])
piy = np.zeros([8000,6000])
Xn = []
Yn = []
Xi = []
Yi = []
num = []
S1,S2,S3,S4,S5,S6,S7,S8,S9,S10=0,0,0,0,0,0,0,0,0,0
T1,T2,T3,T4,T5,T6,T7,T8,T9, T10 = None, None, None, None, None, None, None, None, None, None

for i in range(100000):
  num.append(i)
  if yn == 0:
    yn = 1
  if yi == 0:
    yi = 1
  if xn == 0 and yn == 1:
    T = T1 = -(1/L)*math.log(1-u())
  elif xn == 0:
    T1 = -(1/L)*math.log(1-u())
    T3 = -(1/(H*xi))*math.log(1-u())
    T = min(T1,T3)
  elif yn == 1:
    Pn = (xn*M+yn)/(xn+(yn+yi))
    tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
    taunn = min(C*Pn*xn,M*(N*xn+yn))
    T1 = -(1/L)*math.log(1-u())
    T2 = -(1/(H*xn))*math.log(1-u())
    T4 = -(1/G*yn)*math.log(1-u())
    T9 = -(1/taunn)*math.log(1-u())
    T = min(T1,T2,T4,T9)
  elif xi != 0 and yi < yn:
    tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
    T7 = -(1/tauni)*math.log(1-u())
    T = T7
  elif xi == 0 and yi == 1:
    Pn = (xn*M+yn)/(xn+(yn+yi))
    tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
    taunn = min(C*Pn*xn,M*(N*xn+yn))
    taui = max(C*xi,((xn+xi)*N+(yn+yi)*(xi/(xn+xi))))
    tauC = min(C*K*(1-Pn)*xn,M*K*(M*(xn+xi)+(yn+yi))*(1-Pn))
    T1 = -(1/L)*math.log(1-u())
    T2 = -(1/(H*xn))*math.log(1-u())
    T4 = -(1/G*yn)*math.log(1-u())
    T6 = -(1/(tauC))*math.log(1-u())
    T9 = -(1/taunn)*math.log(1-u())
    T = min(T1,T2,T4,T6,T9)
  elif xi == 0:
    Pn = (xn*M+yn)/(xn+(yn+yi))
    tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
    taunn = min(C*Pn*xn,M*(N*xn+yn))
    taui = max(C*xi,((xn+xi)*N+(yn+yi)*(xi/(xn+xi))))
    tauC = min(C*K*(1-Pn)*xn,M*K*(M*(xn+xi)+(yn+yi))*(1-Pn))
    T1 = -(1/L)*math.log(1-u())
    T2 = -(1/(H*xn))*math.log(1-u())
    T4 = -(1/G*yn)*math.log(1-u())
    T5 = -(1/(G*yi))*math.log(1-u())
    T6 = -(1/(tauC))*math.log(1-u())
    T9 = -(1/taunn)*math.log(1-u())
    T = min(T1,T2,T4,T5,T6,T9)
  else:
    Pn = (xn*M+yn)/(xn+(yn+yi))
    tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
    taunn = min(C*Pn*xn,M*(N*xn+yn))
    taui = min(C*xi,((xn+xi)*N+(yn+yi)*(xi/(xn+xi))))
    tauC = min(C*K*(1-Pn)*xn,M*K*(M*(xn+xi)+(yn+yi))*(1-Pn))
    T1 = -(1/L)*math.log(1-u())
    T2 = -(1/(H*xn))*math.log(1-u())
    T3 = -(1/(H*xi))*math.log(1-u())
    T4 = -(1/G*yn)*math.log(1-u())
    T5 = -(1/(G*yi))*math.log(1-u())
    T6 = -(1/(tauC))*math.log(1-u())
    T7 = -(1/tauni)*math.log(1-u())
    T8 = -(1/taui)*math.log(1-u())
    T9 = -(1/taunn)*math.log(1-u())
    T = min(T1,T2,T3,T4,T5,T6,T7,T8,T9)

  if T == T1:
    pix[xn+1][xi] += T
    piy[yn][yi] += T
    xn += 1
    S1 += 1
  elif T == T2:
    pix[xn-1][xi]
    piy[yn][yi] += T
    xn -= 1
    S2 += 1
  elif T == T3:
    pix[xn][xi-1]
    piy[yn][yi] += T
    xi -= 1
    S3 += 1
  elif T == T4:
    pix[xn][xi] += T
    piy[yn-1][yi] += T
    yn -= 1
    S4 += 1
  elif T == T5:
    pix[xn][xi] += T
    piy[yn][yi-1] += T
    yi -= 1
    S5 += 1
  elif T == T6:
    pix[xn-1][xi+1] += T
    piy[yn][yi] += T
    xn -= 1
    xi += 1
    S6 += 1
  elif T == T7:
    pix[xn-1][xi] += T
    piy[yn][yi+1] += T
    xn -= 1
    yi += 1
    S7 += 1
  elif T == T8:
    pix[xn][xi-1] += T
    piy[yn][yi+1] += T
    xi -= 1
    yi += 1
    S8 += 1
  elif T == T9:
    pix[xn-1][xi] += T
    piy[yn+1][yi] += T
    xn -= 1
    yn += 1
    S9 += 1
  elif T == T10:
    pix[xn][xi] += T
    piy[yn+1][yi] += T
    yn += 1
    S10 += 1
  Xn.append(xn)
  Yn.append(yn)
  Xi.append(xi)
  Yi.append(yi)
  print(xn, yn, xi, yi)
  T1,T2,T3,T4,T5,T6,T7,T8,T9, T10 = None, None, None, None, None, None, None, None, None, None
  Tsim += T
print(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10)
TXn = pd.Series(Xn)
TYn = pd.Series(Yn)
TXi = pd.Series(Xi)
TYi = pd.Series(Yi)

fig, axs = plt.subplots(2, 2, figsize=(10,8))
axs[0, 0].plot(TXn.index,TXn.values)
axs[0, 0].axhline(y=sum(Xn)/100000,c='k')
axs[0, 0].set_title('Leeches sanos')
axs[0, 1].plot(TYn.index,TYn.values, 'tab:orange')
axs[0, 1].axhline(y=sum(Yn)/100000,c='k')
axs[0, 1].set_title('Seeds sanos')
axs[1, 0].plot(TXi.index,TXi.values, 'tab:green')
axs[1, 0].axhline(y=sum(Xi)/100000,c='k')
axs[1, 0].set_title('Leeches infectados')
axs[1, 1].plot(TYi.index,TYi.values, 'tab:red')
axs[1, 1].axhline(y=sum(Yi)/100000,c='k')
axs[1, 1].set_title('Seeds infectados')
fig.tight_layout()
print('Leeches sanos:',sum(Xn)/100000)
print('Seeds sanos:',sum(Yn)/100000)
print('Leeches infectados:',sum(Xi)/100000)
print('Seeds infectados:',sum(Yi)/100000)
print('Valores usados: C ',C,'L',L,'M',M,'H',H,'G',G,'N',N,'K',K)
