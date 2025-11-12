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
h=5
tspan = np.linspace(0, 6, 300)
u0 = [-(np.arccos(1)-np.arccos(h/(L-aL))),0.0]

# %% Define derivative function
def f(t, u, K, use_K):
    phi,w=u
    Ia = M*(aL)**2
    if use_K:
        Ia += M*K**2
    Is=m*(L/2-aL)**2+(m*L**2)/12
    I=Ia+Is
    Ma=M*g*np.sin(phi)*aL
    Ms=m*g*np.sin(phi)*(L/2-aL)
    Mm=Ma-Ms
    dudt = [w,Mm/I]
    return dudt

# %% Solve differential equation for cannot rotate (use K)
sol_cannot = solve_ivp(lambda t, u: f(t, u, K, True), [tspan[0], tspan[-1]], y0=u0, t_eval=tspan, rtol=1e-7)

# %% Solve differential equation for can rotate (no K)
sol_can = solve_ivp(lambda t, u: f(t, u, K, False), [tspan[0], tspan[-1]], y0=u0, t_eval=tspan, rtol=1e-7)

def print_states_at_times(sols, times, labels=None):
    # header
    print("Time (s)\tphi (deg)\tphi_dot (deg/s)")
    for i, tq in enumerate(times):
        sol = sols[i]
        # interpolera till tid tq
        phi = np.interp(tq, sol.t, sol.y[0])
        w   = np.interp(tq, sol.t, sol.y[1])
        print(f"{tq:.2f}\t\t{np.degrees(phi):.3f}\t\t{np.degrees(w):.3f}")
        if labels is not None:
            print(labels[i])

# Exempel: skriv ut de två rader du frågade om
sols = [sol_cannot, sol_can]
print_states_at_times(sols, [2.58, 2.33],
                      ["The mass CANNOT rotate.", "The mass CAN rotate."])

# ...rest of the code if needed...