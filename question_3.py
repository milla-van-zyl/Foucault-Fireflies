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
