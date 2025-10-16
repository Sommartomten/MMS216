import numpy as np
from sympy import symbols, solve, cos, sin, simplify

# declare symbols (using simple names to avoid unicode issues in some environments)
a, N, a_x, a_y, L, M, g, m, m_p, l, phi, phi_1, phi_2, theta, theta_1, theta_2, I = symbols(
    'a N a_x a_y L M g m m_p l phi phi_1 phi_2 theta theta_1 theta_2 I'
)

# Equations as provided (only the first two contain theta_2 and phi_2)
eq1 = (a - 1)*L*(cos(phi)*phi_1**2 + sin(phi)*phi_2) - l*(theta_1 + phi_1)**2*cos(phi + theta) - l*sin(phi + theta)*(theta_2 + phi_2) - a_x
eq2 = (a - 1)*L*(cos(phi)*phi_2 - sin(phi)*phi_1**2) + (theta_1 + phi_1)**2*(-sin(phi + theta)) + cos(phi + theta)*(theta_2 + phi_2) - a_y

# The remaining three equations (they introduce N, a_x, a_y but not theta_2/phi_2)
eq3 = L*m*g*(a*cos(phi)-(1/2-a)) + (1-a)*L*N*sin(theta) - I*phi_1
eq4 = M*(a*L)**2 + m*L*(1/2 - a)**2 + 1/12*m*L**2 - I
eq5 = -m_p*(g*sin(phi+theta) + a_x*cos(phi+theta) + a_y*sin(phi+theta)) - N

# The remaining equations don't contain theta_2 or phi_2, so we only need eq1 and eq2 to solve for phi_2 and theta_2

# Introduce t = theta_2 + phi_2 to simplify linear appearance
t = symbols('t')

# Strategy:
# 1) Solve eq3 and eq5 for N and a_x (or a_y), then substitute into eq1/eq2.
# 2) Introduce t = theta_2 + phi_2 and solve for phi_2 and t.

# Attempt to solve eq3 and eq5 for N and a_x/a_y where possible
sol_N = solve([eq3, eq5], [N, a_x], dict=True)

if sol_N:
    # Substitute N and a_x from solution into eq1 and eq2; also solve eq5 for a_y if needed
    subs_dict = {}
    # pick first solution
    sN = sol_N[0]
    subs_dict.update({N: sN[N], a_x: sN.get(a_x, a_x)})
    # Now try to solve eq5 for a_y explicitly and add to subs
    try:
        sol_ay = solve(eq5.subs(N, sN[N]), a_y, dict=True)
        if sol_ay:
            subs_dict[a_y] = sol_ay[0][a_y]
    except Exception:
        pass

    eq1_sub = eq1.subs(subs_dict)
    eq2_sub = eq2.subs(subs_dict)
else:
    # Fallback: keep a_x,a_y,N symbolic and proceed with elimination later
    eq1_sub = eq1
    eq2_sub = eq2

# Substitute t = theta_2 + phi_2
eq1_t = eq1_sub.subs(theta_2 + phi_2, t)
eq2_t = eq2_sub.subs(theta_2 + phi_2, t)

# Solve for phi_2 and t
sol = solve([eq1_t, eq2_t], [phi_2, t], dict=True)

if not sol:
    print('No symbolic solution found for phi_2 and theta_2 (via t)')
else:
    for s in sol:
        # recover theta_2 = t - phi_2
        s_theta2 = s[t] - s[phi_2]
        print('phi_2 =', simplify(s[phi_2]))
        print('theta_2 =', simplify(s_theta2))