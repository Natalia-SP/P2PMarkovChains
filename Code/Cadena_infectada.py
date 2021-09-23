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
pix = Array{Float64}(undef, 5000,6000)
piy = Array{Float64}(undef, 8000,6000)
Xn = []
Yn = []
Xi = []
Yi = []
num = []
S1,S2,S3,S4,S5,S6,S7,S8,S9,S10=0,0,0,0,0,0,0,0,0,0
T1,T2,T3,T4,T5,T6,T7,T8,T9,T10 = Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing

for i in 1:100000
    append!(num,i)
    if yn == 0
        yn = 1
    end
    if yi == 0
        yi = 1
    end
    if xn == 0 && yn == 1
        T = T1 = -(1/L)*log(1-u())
    elseif xn == 0
        T1 = -(1/L)*log(1-u())
        T3 = -(1/(H*xi))*log(1-u())
        T = min(T1,T3)
    elseif yn == 1
        Pn = (xn*M+yn)/(xn+(yn+yi))
        tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
        taunn = min(C*Pn*xn,M*(N*xn+yn))
        T1 = -(1/L)*log(1-u())
        T2 = -(1/(H*xn))*log(1-u())
        T4 = -(1/G*yn)*log(1-u())
        T9 = -(1/taunn)*log(1-u())
        T = min(T1,T2,T4,T9)
    elseif xi != 0 && yi < yn
        Pn = (xn*M+yn)/(xn+(yn+yi))
        tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
        T7 = -(1/tauni)*log(1-u())
        T = T7
    elseif xi == 0 && yi == 1
        Pn = (xn*M+yn)/(xn+(yn+yi))
        tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
        taunn = min(C*Pn*xn,M*(N*xn+yn))
        taui = max(C*xi,((xn+xi)*N+(yn+yi)*(xi/(xn+xi))))
        tauC = min(C*K*(1-Pn)*xn,M*K*(M*(xn+xi)+(yn+yi))*(1-Pn))
        T1 = -(1/L)*log(1-u())
        T2 = -(1/(H*xn))*log(1-u())
        T4 = -(1/G*yn)*log(1-u())
        T6 = -(1/(tauC))*log(1-u())
        T9 = -(1/taunn)*log(1-u())
        T = min(T1,T2,T4,T6,T9)
    elseif xi == 0
        Pn = (xn*M+yn)/(xn+(yn+yi))
        tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
        taunn = min(C*Pn*xn,M*(N*xn+yn))
        taui = max(C*xi,((xn+xi)*N+(yn+yi)*(xi/(xn+xi))))
        tauC = min(C*K*(1-Pn)*xn,M*K*(M*(xn+xi)+(yn+yi))*(1-Pn))
        T1 = -(1/L)*log(1-u())
        T2 = -(1/(H*xn))*log(1-u())
        T4 = -(1/G*yn)*log(1-u())
        T5 = -(1/(G*yi))*log(1-u())
        T6 = -(1/(tauC))*log(1-u())
        T9 = -(1/taunn)*log(1-u())
        T = min(T1,T2,T4,T5,T6,T9)
    else
        Pn = (xn*M+yn)/(xn+(yn+yi))
        tauni = min(C*(1-Pn)*xi,M*(N*xn+yn)*(1-Pn))
        taunn = min(C*Pn*xn,M*(N*xn+yn))
        taui = min(C*xi,((xn+xi)*N+(yn+yi)*(xi/(xn+xi))))
        tauC = min(C*K*(1-Pn)*xn,M*K*(M*(xn+xi)+(yn+yi))*(1-Pn))
        T1 = -(1/L)*log(1-u())
        T2 = -(1/(H*xn))*log(1-u())
        T3 = -(1/(H*xi))*log(1-u())
        T4 = -(1/G*yn)*log(1-u())
        T5 = -(1/(G*yi))*log(1-u())
        T6 = -(1/(tauC))*log(1-u())
        T7 = -(1/tauni)*log(1-u())
        T8 = -(1/taui)*log(1-u())
        T9 = -(1/taunn)*log(1-u())
        T = min(T1,T2,T3,T4,T5,T6,T7,T8,T9)
    end

    if T == T1
        pix[1+xn+1,1+xi] += T
        piy[1+yn,1+yi] += T
        xn += 1
        S1 += 1
    elseif T == T2
        pix[1+xn-1,1+xi] += T
        piy[1+yn,1+yi] += T
        xn -= 1
        S2 += 1
    elseif T == T3
        pix[1+xn,1+xi-1]
        piy[1+yn,1+yi] += T
        xi -= 1
        S3 += 1
    elseif T == T4
        pix[1+xn,1+xi] += T
        piy[1+yn-1,1+yi] += T
        yn -= 1
        S4 += 1
    elseif T == T5
        pix[1+xn,1+xi] += T
        piy[1+yn,1+yi-1] += T
        yi -= 1
        S5 += 1
    elseif T == T6
        pix[1+xn-1,1+xi+1] += T
        piy[1+yn,1+yi] += T
        xn -= 1
        xi += 1
        S6 += 1
    elseif T == T7
        pix[1+xn-1,1+xi] += T
        piy[1+yn,1+yi+1] += T
        xn -= 1
        yi += 1
        S7 += 1
    elseif T == T8
        pix[1+xn,1+xi-1] += T
        piy[1+yn,1+yi+1] += T
        xi -= 1
        yi += 1
        S8 += 1
    elseif T == T9
        pix[1+xn-1,1+xi] += T
        piy[1+yn+1,1+yi] += T
        xn -= 1
        yn += 1
        S9 += 1
    elseif T == T10
        pix[1+xn,1+xi] += T
        piy[1+yn+1,1+yi] += T
        yn += 1
        S10 += 1
    end
    append!(Xn,xn)
    append!(Yn,yn)
    append!(Xi,xi)
    append!(Yi,yi)
    T1,T2,T3,T4,T5,T6,T7,T8,T9, T10 = Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing
    Tsim += T
end
print(S1,'\t',S2,'\t',S3,'\t',S4,'\t',S5,'\t',S6,'\t',S7,'\t',S8,'\t',S9,'\t',S10)
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
