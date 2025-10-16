from matplotlib import pyplot as plt
import numpy as np
import scipy

#konstanter



def equal_F(x):
    return mj*(1-x)**2-mm*(x**2)
# startläge
r0= 1.5*rj  # utgångs höjd från jordens yta
v0= [9050]  # utgångshastighet ur jordens omloppsbana
vin0=[230]  # start läge på sonden som grader i start på x axeln
temp=[]               
for item in vin0:
    temp.append(np.deg2rad(item))
vin0=temp

def func(t,u):
    x,vx,y,vyZZ=u
    #moon position:
    x_moon = am*np.cos(2*np.pi*t/T)
    y_moon = am*np.sin(2*np.pi*t/T)
    sond_moon= np.sqrt((x_moon-x)**2+(y_moon-y)**2) # avstånd sond till månen

    earth_sond=np.sqrt(x**2+y**2) #distance earth to sond

    dxdv=-G*(mj*x/earth_sond**3 + (mm*(x-x_moon))/ (sond_moon**(3)))

    dydv=-G*(mj*y/earth_sond**3 +(mm*(y-y_moon))/ (sond_moon**(3)))
    return([vx,dxdv,vy,dydv])

# Tidsintervall och utvärderingspunkter
t_span = (0, 10 * 24*3600)  # 10 dagar i sekunder
t_eval = np.linspace(t_span[0], t_span[1], 1000)


def calculate():
# Lösning av ODE-systemet
    for vin in vin0:
        for speed in v0:
            x0=r0*np.cos(vin)
            y0=r0*np.sin(vin)
            vx0=speed*np.cos(vin+np.pi/2)
            vy0=speed*np.sin(vin+np.pi/2)
            sol = scipy.integrate.solve_ivp(func, t_span, [x0, vx0, y0, vy0], t_eval=t_eval)

            x = sol.y[0]
            y = sol.y[2]

            # Plotta banan
            plt.figure(figsize=(8, 8))
            plt.plot(x,y, label=f'sondens bana')
            plt.plot(am * np.cos(2 * np.pi * sol.t / T), am * np.sin(2 * np.pi * sol.t / T), '--', label="Månens bana")
            plt.plot(0, 0, 'yo', label="Jorden")
            plt.xlabel('x-position (m)')
            plt.ylabel('y-position (m)')
            plt.title(f'Bana för objektet över tid med utgång {round(np.rad2deg(vin))} grader och {speed}m/s')
            plt.grid(True)
            plt.axis('equal')
            #plt.savefig(f'bana({round(np.rad2deg(vin))}-{speed}).png')
            xm = am * np.cos(2 * np.pi * sol.t / T)
            ym = am * np.sin(2 * np.pi * sol.t / T)
            distance_to_moon = np.sqrt((xm-sol.y[0])**2 + (ym-sol.y[2])**2)

            plt.figure()
            plt.plot(sol.t / 3600, distance_to_moon / 1e6)
            plt.xlabel("Tid [h]")
            plt.ylabel("Avstånd till månen [Mm]")
            plt.title(f'Avstånd mellan sond och måne över tid med utgång {round(np.rad2deg(vin))} grader och {speed}m/s')
            plt.grid()
        # plt.savefig(f'avstand({round(np.rad2deg(vin))}-{speed}).png')
            plt.show()

print("andel av månens avstånd för netto noll påverkan",scipy.optimize.fsolve(equal_F,0.5))
print(f'omloppshastighet för {r0*10**-3} km över jorden är {round(np.sqrt(G*mj/r0),1)} m/s')
calculate()
print("done")