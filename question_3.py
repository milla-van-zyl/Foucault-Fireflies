from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# Parameters

g = 9.81 # m/s^2
L = 67 # metres
gamma = 0.01 # Friction coefficient

def pendulum(t, state):
	x, y, u_x, u_y = state
	dxdt = u_x
	dydt = u_y
	du_xdt = -(g/L) * x - gamma * u_x
	du_ydt = -(g/L) * y - gamma * u_y

	return [dxdt, dydt, du_xdt, du_ydt]

timespan = (0, 360)

# Initial conditions (x, y, u_x, u_y) which we will change for next questions

initial_state = [1.0, 0.0, 0.0, 0.0]

sol = solve_ivp(pendulum, timespan, initial_state, t_eval = np.linspace(0, 360, 10000))

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

# Plotting x(t) and y(t)

plt.figure(figsize=(10,4))
plt.plot(t, x, label="x(t)")
plt.plot(t, y, label="y(t)")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.title("Pendulum displacement against time")
plt.legend()
plt.grid(True)
plt.show()

# Right now, the initial amplitude is quite large (x = 1), so we only need to change the x initial condition to be smaller and then plot the graphs on latex!

