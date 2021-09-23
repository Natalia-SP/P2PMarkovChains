# Librerías
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def u():
  u = random.uniform(0,1000000)
  return u 

# Estado inicial

i = 0 # Peers Leeches
x = 1  #Peers Seeds

# Parámetros

C = 0.02
L = 1
M = 0.00125
H = 0.01
G = 0.01
N = 0.85
O = 1
P = 1

# Inicializar listas, arrays y algunas variables
TCsim = 0 # Almacenará el tiempo de simulación
pi = Array{Float64}(undef, 200,200) #Almacena el tiempo en que paso por determinado estado
I = [] # Guarda en cada iteración el valor actual de i
X2 = [] # Guarda en cada iteración el valor actual de x
num = [] #Guarda una lista con los números de iteraciones
T1,T2,T3,T4,T5,T6,T7,T8,T9, T10 = Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing

# Ciclo de iteraciones de 1 a 100,000 con un paso de uno

for j in 1:1:100000
    append!(num,j)
    if x == 0
        x = 1
    end
    if i == 0 && x == 1
        T = T1 = -(1/L)*log(1-u())
    end
    if i == 0
        T1 = -(1/L)*log(1-u())
        T3 = -(1/(G*x))*log(1-u())
        T = min(T1,T3)
    elseif x == 1
        tau = min(M*(N*i+x),C*i)
        T1 = -(1/L)*log(1-u())
        T2 = -(1/(H*i))*log(1-u())
        T4 = -(1/tau)*log(1-u())
        T = min(T1,T2,T4)
    else
        tau = min(M*(N*i+x),C*i)
        T1 = -(1/L)*log(1-u())
        T2 = -(1/(H*i))*log(1-u())
        T3 = -(1/(G*x))*log(1-u())
        T4 = -(1/tau)*log(1-u())
        T5 = -(1/P)*log(1-u())
        T6 = -(1/O)*log(1-u())
        T = min(T1,T2,T3,T4,T5,T6)
    end
    
    # Una vez que se obtuvo el tiempo minimo se encuentra
    # al estado que pertenece y conforme a ese estado se
    # modifica ya sea los seeds, leeches o ambos.
    
    if T == T1
        pi[i+1,x] += T
        i += 1
    elseif T == T2
        pi[i-1,x] += T
        i -= 1
    elseif T == T3
        pi[i-1,x-1] += T
        i -= 1
        x -= 1
    elseif T == T4
        pi[i+1,x+1] += T
        i += 1
        x += 1
    elseif T == T5
        pi[i,x] += T
    elseif T == T6
        pi[i+1,x+1] += T
        i += 1
        x += 1
    end
    append!(I,i)
    append!(X2,x)
    println(i, '\t', x)
    TCsim += T
end

# Se grafican los resultados obtenidos
p1 = plot(num,I, label= :none, title="Seeds")
p2 = plot(num,X2, color = :green, label= :none, title="Leeches")
plot(p1, p2, layout = (2, 1))
