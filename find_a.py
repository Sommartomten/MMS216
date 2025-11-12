import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button  # added widgets import
# %% Define initial values, and constants
L = 15.0 #length of arm
g = 9.82 #gravity m/s^2
hight = 4.0 #hight of pole
M = 500.0 #mass of counterweight
m = 50.0 #mass of arm
m_p = 5.0 #mass of projectile
l = 5.0 #length of sling
tspan = np.linspace(0, 2.0, 400)



# sanity print
#print(f"start constraint check: y={ (L - (L*a)) * np.sin(phi0)+ l * np.sin(theta0 + phi0)} ,x={(L - (L*a)) * np.cos(phi0) + l * np.cos(theta0 + phi0)}")




# %% Define derivative function using the extended-model linear system
def f(t, u):


    phi, theta, phidot, thetadot = u

    # For clarity name local variables matching symbolic eqns:

    
    # Kinematic and dynamic equations from step 2
    # Calculate normal force N from eq5
    N = -m_p * (-g * np.sin(phi + theta)    +l*(phidot+thetadot)**2   -L*(1-a)*phidot**2*np.cos(theta))  # Since a_x = a_y = 0 for fixed pivot
    
    # Calculate phi_ddot from eq3 (moment equation)
    moment_motvikt = M * g * np.cos(phi) * (a * L)
    moment_arm = -(m * g * np.cos(phi) * (L/2 - a * L))
    
    
    # Solve for phi_ddot from moment equation
    phidd =(-m*g*(L/2-a*L)*np.cos(phi)+M*g*a*L*np.cos(phi) +(L-a*L)*m_p*np.sin(theta)*(np.cos(theta)*(L-a*L)*phidot**2+l*(phidot+thetadot)**2-g*np.sin(phi+theta)))/(I+m_p*((L-a*L)**2)*np.sin(theta)**2)
    # Calculate theta_ddot using constraint equation
    # Using the simplified form from your original commented code which worked correctly
    thetadd = (((-phidot**2)*(L-a*L)*np.sin(theta)-phidd*(L-a*L)*np.cos(theta)-g*np.cos(phi+theta))/l)-phidd 
    
    dudt = (phidot, thetadot, phidd, thetadd)
    return dudt

# %% Solve differential equation
for a in np.linspace(0.0001,hight/L,100):
    print(f'Testing a={a}, this many percent done: {a/(L/hight)*100}%')
    max=[0,0]
    # Moment of inertia 
    I = ((1/12)*m*(L**2)+m*((L/2)-a*L)**2)+M*(a*L)**2 

    # Initial angles 
    phi0 = (np.arcsin(-hight / (L - (L*a))))  
    theta0 = ( -phi0-np.pi)  # Set theta0 to -pi (180 degrees) for negative x-direction
    phidot0 = 0.0
    thetadot0 = 0.0
    u0 = (phi0, theta0, phidot0, thetadot0)

    sol = solve_ivp(lambda t, u: f(t, u), [tspan[0], tspan[-1]], y0=u0, t_eval=tspan, rtol=1e-7)
    for item in sol.y:
        acceleration=np.sqrt(((a - 1)*L*np.cos(item[0])*item[2] + (item[3] + item[2])*np.cos(item[0] + item[1]))**2+((a - 1)*L*np.sin(item[0])*item[2] + (item[3] + item[2])*np.sin(item[0] + item[1]))**2)

        if acceleration>max[0]:
            max=[acceleration,a]

print(hight/L)
print(f"Max acceleration: {max[0]} at a={max[1]}")