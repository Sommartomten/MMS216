import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
# %% Define initial values, and constants
# lenghts in m, masses in kg
L = 15.0 #length of arm
g = -9.82 #gravity m/s^2
hight = 4.0 #hight of pole
M = 500.0 #mass of counterweight
m = 50.0 #mass of arm
m_p = 5.0 #mass of projectile
l = 5.0 #length of sling
a = 0.25 #relative position of counterweight along arm (0-1)

# Moment of inertia 
I = - M * (L*a)**2 + m * (L/2 - (L*a))**2 + 1/12 * m * L**2

# Initial angles 
phi0 = (np.arcsin(-hight / (L - (L*a))))  
theta0 = (np.pi - phi0)  # Set theta0 to -pi (180 degrees) for negative x-direction

phidot0 = 0.0
thetadot0 = 0.0
u0 = (phi0, theta0, phidot0, thetadot0)

# sanity print
#print(f"start constraint check: y={ (L - (L*a)) * np.sin(phi0)+ l * np.sin(theta0 + phi0)} ,x={(L - (L*a)) * np.cos(phi0) + l * np.cos(theta0 + phi0)}")

tspan = np.linspace(0, 2.0, 400)

# %% Define derivative function using the extended-model linear system
def f(t, u):


    phi, theta, phidot, thetadot = u

    # For clarity name local variables matching symbolic eqns:

    
    # Kinematic and dynamic equations from step 2
    # Calculate normal force N from eq5
    N = -m_p * (g * np.sin(phi + theta)-l*(phidot+thetadot)**2-L*(1-a)*phidot**2*np.cos(theta))  # Since a_x = a_y = 0 for fixed pivot
    
    # Calculate phi_ddot from eq3 (moment equation)
    moment_motvikt = M * g * np.cos(phi) * (a * L)
    moment_arm = -(m * g * np.cos(phi) * (L/2 - a * L))
    moment_inertia = I
    
    # Solve for phi_ddot from moment equation
    phidd = (moment_motvikt + moment_arm + (1-a)*L*N*np.sin(theta)) / (moment_inertia+m_p*(L-L*a)**2*np.sin(theta)**2)
    
    # Calculate theta_ddot using constraint equation
    # Using the simplified form from your original commented code which worked correctly
    thetadd = ((-phidot**2*L*(1-a)*np.sin(theta) - phidd*L*(1-a)*np.cos(theta)-g*np.cos(phi+theta))/l)-phidd
    
    dudt = (phidot, thetadot, phidd, thetadd)
    return dudt

# %% Solve differential equation
sol = solve_ivp(lambda t, u: f(t, u), [tspan[0], tspan[-1]], y0=u0, t_eval=tspan, rtol=1e-7)

## Extract data
time = sol.t
phi = np.degrees(sol.y[0])     # converting from radians to degrees for plotting
theta = np.degrees(sol.y[1])

# Angular accelerations/velocities: Note solve_ivp stores derivatives as states we returned
omega_phi = np.degrees(sol.y[2])   # phidot in degrees/s (state index 2 is phi_dot)
omega_theta = np.degrees(sol.y[3]) # theta_dot in degrees/s

# Calculate trajectory assuming arm and sling lengths (use radians for trig)
phi_rad = sol.y[0]
theta_rad = sol.y[1]

X = (L - (L*a)) * np.cos(phi_rad) + l * np.cos(theta_rad + phi_rad)
Y = (L - (L*a)) * np.sin(phi_rad) + l * np.sin(theta_rad + phi_rad)

x_arm = (L - (L*a)) * np.cos(phi_rad)
y_arm = (L - (L*a)) * np.sin(phi_rad)
x_sling = x_arm + l * np.cos(theta_rad + phi_rad)
y_sling = y_arm + l * np.sin(theta_rad + phi_rad)

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

axs[0, 0].plot(time, phi)
axs[0, 0].set_title('Rotating angle of the arm')
axs[0, 0].set_xlabel('Time (s)')
axs[0, 0].set_ylabel('phi (Degrees)')

axs[0, 1].plot(time, theta, color='orange')
axs[0, 1].set_title('Rotating angle of the sling')
axs[0, 1].set_xlabel('Time (s)')
axs[0, 1].set_ylabel('theta (Degrees)')

axs[1, 0].plot(time, omega_phi, color='green')
axs[1, 0].set_title('Angular velocity of the arm')
axs[1, 0].set_xlabel('Time (s)')
axs[1, 0].set_ylabel('omega_phi (Degrees/s)')

axs[1, 1].plot(time, omega_theta, color='red')
axs[1, 1].set_title('Angular velocity of the sling')
axs[1, 1].set_xlabel('Time (s)')
axs[1, 1].set_ylabel('omega_theta (Degrees/s)')

plt.tight_layout()

# Trajectory plot
plt.figure(figsize=(6, 6))
plt.plot(X, Y)
plt.scatter(X[0], Y[0], color='black', label='Start')
plt.scatter(X[-1], Y[-1], color='yellow', label='End')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Trajectory (zero point is at the bottom of the support)')
plt.legend()
plt.show()

# Setup figure and axis for animation
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-((L - (L*a)) + l + 2), (L - (L*a)) + l + 2)
ax.set_ylim(-((L - (L*a)) + l + 2), (L - (L*a)) + l + 2)
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2, markersize=8)

def init():
    line.set_data([], [])
    return line,

def update(frame):
    # frame corresponds to time index
    x = [0, x_arm[frame], x_sling[frame]]
    y = [0, y_arm[frame], y_sling[frame]]
    line.set_data(x, y)
    return line,

# Create animation
ani = FuncAnimation(fig, update, frames=len(tspan), init_func=init, blit=True, interval=20)

plt.title("Double Pendulum Animation")
plt.show()