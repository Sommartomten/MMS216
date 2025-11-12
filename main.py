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
<<<<<<< HEAD
tspan = np.linspace(0, 6, 300)
u0 = [-(np.arccos(1)-np.arccos(h/(L-aL))),0.0]
=======
tspan = np.linspace(0, 2.58, 300)
u0 = [np.arcsin(((-h)/(L-aL))),0.0]
>>>>>>> f06a9a1a73bfca92f7b51b72a5035f3ddbbbc3ab


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

def print_states_at_times(sol, times, labels=None):
    # header
    print("Time (s)\tphi (deg)\tphi_dot (deg/s)")
    for i, tq in enumerate(times):
        # interpolera till tid tq
        phi = np.interp(tq, sol.t, sol.y[0])
        w   = np.interp(tq, sol.t, sol.y[1])
        print(f"{tq:.2f}\t\t{np.degrees(phi):.3f}\t\t{np.degrees(w):.3f}")
        if labels is not None:
            print(labels[i])

# Exempel: skriv ut de två rader du frågade om
print_states_at_times(sol, [2.58, 2.33],
                      ["The mass CANNOT rotate.", "The mass CAN rotate."])
# ...existing code...
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

# --- Tillagt: energiberäkningar (T, V, E=T+V) ---
# beräkna tröghetsmoment (samma som i f)
Ia = M*(aL)**2 + M*K**2
Is = m*(L/2 - aL)**2 + (m*L**2)/12
I = Ia + Is

phi = sol.y[0, :]
w = sol.y[1, :]

T = 0.5 * I * w**2                        # kinetisk energi
C = M*g*aL - m*g*(L/2 - aL)               # koefficient för cos(phi) i potentiell energi
V = C * np.cos(phi)                       # potentiell energi (val av referens)
E = T + V                                 # total energi

# Skriv ut avvikelse från första till sista tidpunkten
E0 = E[0]
dE_max = np.max(np.abs(E - E0))
print(f"Max avvikelse i energi från startvärde: {dE_max:.3e} (bör vara nära 0 om konservativt)")

plt.figure(5)
plt.plot(sol.t, T, label='T (kinetisk)')
plt.plot(sol.t, V, label='V (potentiell)')
plt.plot(sol.t, E, label='E (total)')
plt.legend()
plt.xlabel('t')
plt.ylabel('Energi')
plt.title('Energi: T, V och total (kontroll av T1+V1 = T2+V2)')
plt.show(block=True)