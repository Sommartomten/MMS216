import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
# %% Define time spans, initial values, and constants
#lenghts in M, masses in KG
length=15
g=9.82
hight=4
Mass_motvikt=500
Mass_length=50
Mass_partikel=5
Length_string=5
a=0.25
length_counterwight=length*a

I=Mass_motvikt*length_counterwight**2+Mass_length*(length/2-length_counterwight)**2+1/12*Mass_length*length**2
phi0=-(np.arcsin(hight/(length-length_counterwight))-np.arccos(1))
theta0=-((-phi0)-np.arccos(-1))


phidot0=0
thetadot0=0
u0 = (phi0,theta0,0,0)
print((length-length_counterwight) * np.sin(phi0)+ Length_string * np.sin(theta0+phi0))

tspan = np.linspace(0, 2.0, 500)


#a,N,a_x,a_y,L,M,g,m,m_p,l,ϕ,ϕ_1,ϕ_2,θ,θ_1,θ_2,I
#a,N,B,   C, L,M,g,m,D,  l,Z, E,  F, W, X,  Y,I
# %% Define derivative function
def f(t, u):
    phi, theta,phidot,thetadot=u
    #x=(length-length_counterwight)*np.sin(phi)+Length_string*np.sin(phi+theta)
    #y=(length-length_counterwight)*np.cos(phi)+Length_string*np.cos(phi+theta)
    #phidd=(np.sin(2*phi+theta)*phidot**2 + (Length_string/(length-length*a))*(phidot+thetadot)**2 + g*np.sin(phi+theta)/(length-length*a)) / np.cos(2*phi+theta)
    #thetadd=np.sin(2*phi+theta)*thetadot**2 #(-length*a*Mass_motvikt*g*np.sin(phi) +Mass_length*g*length/2 -Mass_length*g*length*a + Mass_partikel*g*(length-length*a)*np.cos(theta)*np.sin(phi+theta)-Mass_partikel*(length-length*a)**2*np.cos(theta)*np.sin(2*phi+theta)*phidot**2-2*Mass_partikel*Length_string*(length-length*a)*np.cos(theta)*np.sin(phi+theta)*np.cos(phi+theta)*(phidot+thetadot)**2-Mass_partikel*(length-length*a)*np.cos(theta)*((length-length*a)*np.cos(2*phi+theta) + Length_string*(np.cos(phi+theta) + np.sin(phi+theta)))*phidd )/(Mass_partikel*Length_string*(length-length*a)*np.cos(theta)*(np.cos(phi+theta) + np.sin(phi+theta)))
    phidd=(Mass_motvikt+(length-length_counterwight)*np.sin(theta)*Mass_partikel*(np.cos(theta)*(length-length_counterwight)*phidot**2+(phidot+thetadot)**2*Length_string-g*np.sin(phi+theta)))/(I+Mass_partikel*(length-length_counterwight)**2*np.sin(theta)*np.sin(phi+theta))
    thetadd=(-phidot**2*(length-length_counterwight)*np.sin(theta)-phidd*(length-length_counterwight)*np.cos(theta))/Length_string

    dudt = (phidot,thetadot,phidd,thetadd)
    return dudt



# %% Solve differential equation
sol = solve_ivp(lambda t, u: f(t, u), [tspan[0], tspan[-1]], y0=u0, t_eval=tspan, rtol=1e-7)


## Extract data
time = sol.t
phi = np.degrees(sol.y[0])     # converting from radians to degrees if needed
theta = np.degrees(sol.y[1])

# Angular velocities from solve_ivp solution (already derivatives)
omega_phi = np.degrees(sol.y[2])  # phi_dot in degrees/s
omega_theta = np.degrees(sol.y[3])  # theta_dot in degrees/s

# Calculate trajectory assuming arm and sling lengths

X = (length-length_counterwight) * np.cos(np.radians(phi)) + Length_string * np.cos(np.radians(theta+phi))
Y = (length-length_counterwight) * np.sin(np.radians(phi)) + Length_string *np.sin(np.radians(theta+phi))

x_arm =  (length-length_counterwight) * np.cos(phi)
y_arm =  (length-length_counterwight) * np.sin(phi)
x_sling = x_arm + Length_string * np.cos(theta+phi)
y_sling = y_arm + Length_string * np.sin(theta+phi)

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
ax.set_xlim(-((length-length_counterwight) + Length_string + 2), (length-length_counterwight) + Length_string+ 2)
ax.set_ylim(-((length-length_counterwight) + Length_string + 2), (length-length_counterwight) + Length_string + 2)
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
