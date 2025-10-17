import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button  # added widgets import
# %% Define initial values, and constants
# lenghts in m, masses in kg
L = 15.0 #length of arm
g = 9.82 #gravity m/s^2
hight = 4.0 #hight of pole
M = 500.0 #mass of counterweight
m = 50.0 #mass of arm
m_p = 5.0 #mass of projectile
l = 5.0 #length of sling
a = 0.25 #relative position of counterweight along arm (0-1)

# Moment of inertia 
I = ((1/12)*m*(L**2)+m*((L/2)-a*L)**2)+M*(a*L)**2 

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

# Add slider to control frame index
dt = tspan[1] - tspan[0]  # added: time step in seconds
slider_ax = fig.add_axes([0.15, 0.02, 0.65, 0.03])  # x, y, width, height
# use time (seconds) for slider range and step
frame_slider = Slider(slider_ax, 'Time (s)', tspan[0], tspan[-1], valinit=tspan[0], valstep=dt, valfmt='%0.3f')

# Interaction state for slider dragging
slider_active = False
slider_idx = 0

# Callback to update plot when slider is moved (do NOT stop animation)
def on_slider(val):
    global slider_active, slider_idx
    # val is time in seconds -> convert to nearest frame index
    slider_active = True
    idx = int(round((val - tspan[0]) / dt))
    idx = max(0, min(len(tspan)-1, idx))
    slider_idx = idx
    # Update line data to the selected frame immediately
    x = [0, x_arm[slider_idx], x_sling[slider_idx]]
    y = [0, y_arm[slider_idx], y_sling[slider_idx]]
    line.set_data(x, y)
    fig.canvas.draw_idle()

frame_slider.on_changed(on_slider)

# Detect mouse release to stop "dragging" mode
def on_slider_release(event):
    global slider_active
    if event.inaxes == slider_ax:
        slider_active = False

fig.canvas.mpl_connect('button_release_event', on_slider_release)

# Optional: add Play/Pause button to control animation and keep slider synced
play_ax = fig.add_axes([0.82, 0.02, 0.12, 0.04])
play_button = Button(play_ax, 'Play')

anim_running = True
def toggle_play(event):
    global anim_running
    if anim_running:
        ani.event_source.stop()
        play_button.label.set_text('Play')
    else:
        ani.event_source.start()
        play_button.label.set_text('Pause')
    anim_running = not anim_running

play_button.on_clicked(toggle_play)

# Keep slider synced while animating by updating slider in update()
_orig_update = update
def update_with_slider(frame):
    # If user is dragging the slider, show the slider frame instead of the animation frame
    if slider_active:
        idx = slider_idx
    else:
        idx = int(frame)  # use current animation frame
    result = _orig_update(idx)
    # update slider value only when not interacting (avoid triggering on_slider)
    if not slider_active:
        frame_slider.eventson = False
        frame_slider.set_val(tspan[idx])   # set slider in seconds
        frame_slider.eventson = True
    return result

# swap animation to use the wrapped updater
ani.event_source.stop()
ani._func = update_with_slider
ani.event_source.start()

plt.title("Double Pendulum Animation")
plt.show()