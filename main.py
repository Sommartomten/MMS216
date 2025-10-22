import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# %% Define time spans, initial values, and constants
L=15
aL=0.25*15
g=9.82
K=2
m=50
M=500
h=4
tspan = np.linspace(0, 2.58, 300)
u0 = [np.arcsin(((-h)/(L-aL))),0.0]


# %% Define derivative function

def f(t, u, K):
    phi,w=u
    Ia=(M*(aL)**2+M*K**2)
    Is=m*(L/2-aL)**2+(m*L**2)/12
    I=Ia+Is
    Ma=M*g*np.cos(phi)*aL
    Ms=m*g*np.cos(phi)*(L/2-aL)
    Mm=Ma-Ms
    dudt = [w,Mm/I]
    return dudt
# %% Solve differential equation
sol = solve_ivp(lambda t, u: f(t, u, K), [tspan[0], tspan[-1]], y0=u0, t_eval=tspan, rtol=1e-7)
print(sol)
# %% Plot states
plt.figure(1)
for i in range(sol.y.shape[0]):
    plt.plot(sol.t, sol.y[i])
#plt.show(block=True)
plt.legend(['phi', 'w'])
plt.xlabel('t')
plt.ylabel('tillstånd')
plt.title('lösning som funktion av tid')

# vi kan bestämma Fn direkt från vår tillståndsvektor
Fn = m*L/2*sol.y[1]**2-m*g*np.cos(sol.y[0])

# För Fs är det inte lika enkelt, eftersom den behöver derivata av tillstånden
# vi anropar vår derivatafunktion för det (i en loop)
Fs = 0*Fn
Fn2 = 0*Fn
#for i in range(sol.t):
#print(sol.t)
#print(sol.y)
#for i in range(0,20):
for i in range(sol.y.shape[1]):
    #print(sol.t[i])
    #print(sol.y[0, i])
    #print(i)
    yy = f(sol.t[i], sol.y[:, i], K)

    Fs[i] = m*L/2*yy[1]-m*g*np.sin(sol.y[0, i])
    Fn2[i] = m*L/2*yy[0]**2-m*g*np.cos(sol.y[0, i])
#print(Fn-Fn2)

plt.figure(2)
plt.plot(sol.t, Fs)
plt.plot(sol.t, Fn)
plt.legend(['Fs', 'Fn'])
plt.xlabel('t')
plt.ylabel('F')
#plt.show(block=True)
plt.title('Stödkrafter (Fs,Fn) som funktion av tid')

Fx = Fs*np.cos(sol.y[0])-Fn*np.sin(sol.y[0])
Fy = -Fs*np.sin(sol.y[0])-Fn*np.cos(sol.y[0])


plt.figure(3)
plt.plot(sol.t,Fx)
plt.plot(sol.t,Fy)
plt.legend(['Fx', 'Fy'])
plt.xlabel('t')
plt.ylabel('F')
plt.title('Stödkrafter (Fx,Fy) som funktion av tid')


plt.figure(4)
plt.plot(sol.y[0],Fx)
plt.plot(sol.y[0],Fy)
plt.legend(['Fx', 'Fy'])
plt.xlabel('phi')
plt.ylabel('F')
plt.title('Stödkrafter (Fx,Fy) som funktion av vinkel')
plt.show(block=True)