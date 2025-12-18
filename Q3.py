# from scipy.signal import hilbert ## Needed for the Hilbert Transform Envelope
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# Parameters
m = 28 # kg
g = 9.81 # m/s^2
L = 67 # metres
gamma = 0.01 # Friction coefficient

def pendulum(t, state):
	x, y, u_x, u_y = state
	dxdt = u_x
	dydt = u_y
	du_xdt = -x/L**2 * (L**2*(u_x**2 + u_y**2) - (x*u_y - y*u_x)**2/(L**2 - x**2 - y**2)) - (g/L**2)*x*np.sqrt(L**2 - x**2 - y**2) - (gamma/m)*u_x
	du_ydt = -y/L**2 * (L**2*(u_x**2 + u_y**2) - (x*u_y - y*u_x)**2/(L**2 - x**2 - y**2)) - (g/L**2)*y*np.sqrt(L**2 - x**2 - y**2) - (gamma/m)*u_y

	return [dxdt, dydt, du_xdt, du_ydt]

timespan = (0, 120)

# Initial conditions (x, y, u_x, u_y) which we will change for next questions

initial_state = [1.0, 0.0, 0.0, 0.0]

sol = solve_ivp(pendulum, timespan, initial_state, t_eval = np.linspace(0, 120, 10000))

t = sol.t
x = sol.y[0]
y = sol.y[1]


# Period calculator

def period_estimate(t, x):
	crossings = [] ### Makes an empty list storing times when the pendulum crosses x=0 upward (i.e going from -ve to +ve)
	for i in range(len(x)-1): ### Loops through every point besides the last one as we are looking ahead to x[i+1]
		if x[i] < 0 and x[i+1] >= 0: ### Pendulum is on the -ve side AND the next point is zero or +ve, so must have gone through centre between t[i] & t[i+1]
			crossings.append(t[i])

	if len(crossings) < 2: ### We cannot mesaure a period if there are more than 2 crossings
		return None

	return 2 * np.mean(np.diff(crossings)) ### 1 full period is 2x the "crossings", so we take the mean time between the crossings and x2

print("Estimated period:", period_estimate(t, x))

# Envelope function
#envelope = np.abs(hilbert(x))

# Plotting x(t) and y(t)
plt.figure(figsize=(10,5))
plt.plot(t, x, 'b-', linewidth=3, label='x(t)')
plt.plot(t, y, 'r-', label='y(t)')
#plt.plot(t, envelope, linewidth=2, label="Envelope Function")
plt.xlabel("Time, t (s)", fontsize=30)
plt.ylabel("Displacement, x(t) (m)", fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=30)
plt.grid(True, alpha=0.3)
plt.legend(loc=1, prop={'size': 30})
plt.tight_layout(rect=[0.02, 0, 0.98, 1], pad=0.5)
plt.show()



