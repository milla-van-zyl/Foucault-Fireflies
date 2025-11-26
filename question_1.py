from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# Parameters

g = 9.81 # m/s^2
L = 67 # metres
sigma = np.radians(48.1014) # Converts degrees to radians with numpy
omega = 2 * np.pi / 86400 # Seconds

# The full system, but broken down into 4 first order equations:
# u_x = dx/dt,
# u_y = dy/dt,
# du_x/dt = 2*omega*cos(sigma)*u_y - (g*x)/L,
# du_y/dt = -2*omega*cos(sigma)*u_x - (g*y)/L

def pendulum(t, state): # State of system at some time t
	x, y, u_x, u_y = state
	dxdt = u_x
	dydt = u_y
	du_xdt = 2 * omega * np.cos(sigma) * u_y - (g/L) * x
	du_ydt = -2 * omega * np.cos(sigma) * u_x - (g/L) * y

	return [dxdt, dydt, du_xdt, du_ydt]


timespan = (0, 10 * 3600)

# Initial conditions (x, y, u_x, u_y)

initial_state = [1.0, 0.0, 0.0, 0.0]

# Solving

sol = solve_ivp(pendulum, timespan, initial_state, t_eval = np.linspace(0, 10*3600, 1000000))

x = sol.y[0]
y = sol.y[1]
print("x(t) at first 5 time points:", sol.y[0][:5])
print("y(t) at first 5 time points:", sol.y[1][:5])







