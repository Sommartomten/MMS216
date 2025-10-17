import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy.integrate import solve_ivp

sin = np.sin
cos = np.cos

# Parametrar
motvikt = 0   # 1 innebär fast inspänd vikt, sätt till 0 om rotera kring sin egen axel
L = 15        # Stång längd
a = 0.25      # Motvikt position ratio
h = 4         # A höjd
g = 9.82      # 
m = 50        # Massa för stång
mp = 5        # Partikel massa
l = 5         # Plungans längd
M = 500       # Motvikt massa
k = 2         # Tröghetsradie

#Tröghet kring A
I_stång = ((1/12)*m*(L**2)+m*((L/2)-a*L)**2)
I_motvikt = M*(k**2)*motvikt+M*(a*L)**2 
I_tot = I_stång + I_motvikt


def rörelseekvationer(t, u):
    """
    Rörelseekvationer
    y[0] = φ (angle)
    y[1] = θ (angle)
    y[2] = φ_dot (angular velocity)
    y[3] = θ_dot (angular velocity)
    """
    phi,theta,phi_dot,theta_dot = u

    phi_ddot = (-m*g*(L/2-a*L)*cos(phi)+M*g*a*L*cos(phi) +(L-a*L)*mp*sin(theta)*(cos(theta)*(L-a*L)*phi_dot**2+l*(phi_dot+theta_dot)**2-g*sin(phi+theta)))/(I_tot+mp*((L-a*L)**2)*sin(theta)**2)

    theta_ddot = (((-phi_dot**2)*(L-a*L)*sin(theta)-phi_ddot*(L-a*L)*cos(theta)-g*cos(phi+theta))/l)-phi_ddot 
    
    return [phi_dot, theta_dot, phi_ddot, theta_ddot]


# Begynnelsevillkor
if h>(L-a*L+l):                             #L-a+l är kortare än h, dvs änden på stång nuddar ej marken.
    phi_0 = -np.pi/2
    theta_0 = 0.01
elif h>(L-a*L):                             #L-a är kortare än h, dvs änden på stång nuddar ej marken. Theta negativ
    phi_0 = -np.pi/2                          
    theta_0 = -np.arccos((h-(L-a*L))/l)
elif h == 0:
    phi_0 = 0.01                            #Begynnelsevinkel for armen
    theta_0 = -(np.pi+phi_0)                #Begynnelsevinkel för slungan
else:
    phi_0 = -np.arcsin(h/(L-a*L))           #Begynnelsevinkel for armen
    theta_0 = -(np.pi+phi_0)                #Begynnelsevinkel för slungan

phi_dot_0 = 0.0                             #BEgynnelsevinkalhastighet för armen
theta_dot_0 = 0.0                           #Begynnelsevinkelhastighet för slungan

u0 = [phi_0, theta_0, phi_dot_0, theta_dot_0]

# Tidsspann:
t_end = 2
tspan = np.linspace(0, t_end, 100)

# Differentialekvationslösare:
sol = solve_ivp(rörelseekvationer, [0, t_end], u0, t_eval=tspan, method='RK45', rtol=1e-8)

phi = sol.y[0]
theta = sol.y[1]
phi_dot = sol.y[2]
theta_dot = sol.y[3]

#PLOT 1 Vinklar och vinkelhastigheter ###################################

plt.figure(1)
labels = ['φ(t)', 'θ(t)', 'φ_dot(t)', 'θ_dot(t)']
for i in range(sol.y.shape[0]):
    plt.plot(sol.t, sol.y[i], label=labels[i])
    # Add endpoint annotation
    end_t = sol.t[-1]
    end_val = sol.y[i][-1]
    plt.plot(end_t, end_val, 'o', markersize=5)
    plt.annotate(f'({end_t:.2f}, {(end_val*(180/np.pi)):.3f})', 
                xy=(end_t-0.1, end_val), 
                xytext=(5, 5), 
                textcoords='offset points',
                fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
plt.legend()
plt.xlabel('t')
plt.ylabel('φ, θ')
plt.title('')

# set y ticks in multiples of pi
ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(base=np.pi/8))

def pi_formatter(val, pos):
    if val == 0:
        return '0'
    coeff = val / np.pi
    return f'{coeff:.2f}π'

ax.yaxis.set_major_formatter(FuncFormatter(pi_formatter))

# Create second y-axis for grader
ax2 = ax.twinx()
ax2.set_ylabel('Degrees (°)')
ax2.set_ylim(np.array(ax.get_ylim()) * 180 / np.pi)
ax2.yaxis.set_major_locator(MultipleLocator(base=45))  # Every 45 degrees
ax2.grid(True, alpha=0.2, linestyle='--')

# Create third y-axis för radianer
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))  # Offset the third axis
ax3.set_ylabel('Radians')
ax3.set_ylim(ax.get_ylim())
ax3.yaxis.set_major_locator(MultipleLocator(base=1))  #radian
ax3.grid(True, alpha=0.2, linestyle=':')

max_φ_dot = np.max(sol.y[2])  # Maxvärden
max_θ_dot = np.max(sol.y[3]) 

print(f"Maximal vinkelhastighet φ_dot(t): {max_φ_dot}")
print(f"Maximal vinkelhastighet θ_dot(t): {max_θ_dot}")

###############################################################

#PLOT 2: Rörlig kastmaskin #####################################

# Calculate positions
Ax, Ay = 0, h

# Position of motvikt B (at avstånd a*L from A along stången)
Bx = Ax + -a*L*np.cos(phi)
By = Ay + -a*L*np.sin(phi)

# Position av stång ände C
Cx = Ax + (L-(a*L))*np.cos(phi)
Cy = Ay + (L-(a*L))*np.sin(phi)

# Position av partikel Q
Qx = Ax + l*np.cos(phi+theta)+(L-a*L)*np.cos(phi)
Qy = Ay + l*np.sin(phi+theta)+(L-a*L)*np.sin(phi)

# Set up plot
margin = 0.5*L
xmin = min(Bx.min(), Cx.min()) - margin
xmax = max(Bx.max(), Cx.max()) + margin
ymin = min(By.min(), Cy.min(), Ay) - margin
ymax = max(By.max(), Cy.max(), Ay) + margin

fig, (ax1) = plt.subplots(1, figsize=(5, 5))

#subplot
ax1.set_aspect('equal')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)
ax1.set_title("Swinging Rod with Block")
ax1.set_xlabel("x (m)")
ax1.set_ylabel("y (m)")
ax1.grid(True, alpha=0.3)

#ground
support_width = 0.1
support_x = [-support_width/2, support_width/2, support_width/2, -support_width/2, -support_width/2]
support_y = [0, 0, Ay, Ay, 0]
ax1.plot(support_x, support_y, 'k-', linewidth=2)

#elements
stång, = ax1.plot([], [], 'b-', linewidth=3, label='Stång')
motvikt = plt.Rectangle((0, 0), 1, 1, color='red', label='Motvikt')
ax1.add_patch(motvikt)
vridpunkt, = ax1.plot(Ax, Ay, 'ko', markersize=5, label='Partikel')
slunga, = ax1.plot([], [], 'k-', linewidth=1, label='Slunga')
partikel = plt.Circle((0, 0), radius=0.3, color='green')
ax1.add_patch(partikel)

ax1.legend(loc='upper right')

def animate(frame):
    """Animation function"""
    # Current angle
    current_φ = phi[frame]

    # Update motvikt (as a small rectangle)
    motvikt_size = 1
    motvikt_x = Bx[frame] - motvikt_size/2
    motvikt_y = By[frame] - motvikt_size/2
    motvikt.set_xy((motvikt_x, motvikt_y))

    # Update partikel (as a small boll)
    partikel_x = Qx[frame]
    partikel_y = Qy[frame]
    partikel.set_center((partikel_x, partikel_y))

    # Update stång
    stång_x = [Bx[frame], Cx[frame]]
    stång_y = [By[frame], Cy[frame]]
    stång.set_data(stång_x, stång_y)

    # Update slunga
    slunga_x = [Cx[frame], Qx[frame]]
    slunga_y = [Cy[frame], Qy[frame]]
    slunga.set_data(slunga_x, slunga_y)
    
    # Update title with current values
    ax1.set_title(f"Kastmaskin (t={sol.t[frame]:.2f}s, φ={current_φ/np.pi:.2f} π rad)")
    
    return stång, motvikt, slunga, partikel

anim = FuncAnimation(fig, animate, frames=len(tspan),
                    interval=t_end, blit=False, repeat=True)

# ############################################

# PLOT 3 ################### Plot partikelns bana i xy-planet
plt.figure(figsize=(8, 6))
plt.plot(Qx, Qy, 'b-', linewidth=2, label='Partikelns bana')
plt.plot(Qx[0], Qy[0], 'go', markersize=3, label='Startposition')
plt.plot(Qx[-1], Qy[-1], 'ro', markersize=3, label='Slutposition')
plt.text(Qx[-1], Qy[-1], f'  ({Qx[-1]:.2f}, {Qy[-1]:.2f})', 
         fontsize=10, verticalalignment='bottom')
plt.axhline(y=0, color='black', linestyle='-', linewidth=2, label='')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Partikelns rörelse i xy-planet')
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.legend()

##############################################


plt.show()