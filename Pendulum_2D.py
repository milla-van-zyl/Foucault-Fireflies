#
"""
A sample program to solve the system of equations

dx/dt  = vx
dvx/dt = -goL*x-gamma*vx
dy/dt  = vy
dvy/dt = -goL*y-gamma*vy

@author: Bernard
"""

import numpy as np
import ODE_RK4 
import matplotlib.pyplot as plt

# Some
goL = 9.81/67   # g/L where L is the length of the pendulum
gamma = 0.0 # friction parameter
V0 = [1,0,0.5,1] # the initial condition [x0,vx0,y0,vy0]
dt = 0.001 # time step
fig_dt = 0.01 # plotting interval 
t0 = 0 # start ime
t_end = 60 # end time

#####################################################
# The function describing the right and side of the equation
# dV/dt = F(t,V)
# where V is a vector containing the current position (x,y) and corresponding
# speeds (vx,vy)of the pendulum in the order (x, vx, y, vy)
# It returns the right hand side of teh equation
####################################################
def F(t,V):
    eq = np.array([0,0,0,0], dtype=float)
    x  = V[0]
    vx = V[1]
    y  = V[2]
    vy = V[3]
    eq[0] = vx
    eq[1] = -goL*x-gamma*vx
    eq[2] = vy
    eq[3] = -goL*y-gamma*vy
    
    return eq

if __name__ == "__main__":

  # This is the line that performs the integration of the equation
  # F : the right hand side of the equation
  # V0: the intial condition (x0,vx0,y0,vy0)
  # t0 : the initial time, usuaully 0
  # t_end : when to stop the integration
  # dt : the integration time step.
  # fig_dt : the time intervals between data output for figures.
  # The values returned:
  # t_list : the time values at which data for a graph has been saved
  # V_list : the values saved for graphics: as a list off [x, vx, y, vy ]   
    
  t_list, V_list =  ODE_RK4.runge_kutta_2nd_order_system(F, V0, t0, t_end, dt, fig_dt)


  
  # The code below gerenartes a number of graphs

  x_list = [val[0] for val in V_list]
  min_list, max_list = ODE_RK4.get_min_max(t_list, x_list)
  t_list_min = [val[0] for val in min_list]
  v_list_min = [val[1] for val in min_list]
  t_list_max = [val[0] for val in max_list]
  v_list_max = [val[1] for val in max_list]
  # Plot x[t]
  plt.figure(figsize=(10,4))
  plt.plot(t_list, x_list, 'b-', linewidth=3, label='x(t)')
  plt.plot(t_list_min, v_list_min, 'r-', linewidth=3, label='Local Minima', alpha=0.7)
  plt.plot(t_list_max, v_list_max, 'g-', linewidth=3, label = 'Local Maxima', alpha=0.7)
  plt.xlabel("Time, t (s)", fontsize=18.5)
  plt.ylabel("Displacement, x(t) (m)", fontsize=18.5)
  plt.tick_params(axis='both', which='major', labelsize=15)
  plt.grid(True, alpha=0.3)
  plt.legend(loc=1, prop={'size': 16})
  plt.tight_layout(rect=[0.02, 0, 0.98, 1], pad=0.5)
  plt.show()

#   # Plot min val of x
#   min_list, max_list = ODE_RK4.get_min_max(t_list, x_list)
#   t_list_min = [val[0] for val in min_list]
#   v_list_min = [val[1] for val in min_list]
#   plt.plot(t_list_min, v_list_min)
#   plt.xlabel("t", fontsize=20)
#   plt.ylabel("min(x(t))", fontsize=20)
#   plt.tight_layout(rect=[0.02, 0, 0.98, 1], pad=0.5)
#   plt.show()
# 
#   # Plot max val of x
#   t_list_max = [val[0] for val in max_list]
#   v_list_max = [val[1] for val in max_list]
#   plt.plot(t_list_max, v_list_max)
#   plt.xlabel("t", fontsize=20)
#   plt.ylabel("max(x(t))", fontsize=20)
#   plt.tight_layout(rect=[0.02, 0, 0.98, 1], pad=0.5)
#   plt.show()

#   # Plot min and max on the same figure
#   plt.figure(figsize=(10,6))
#   plt.plot(t_list_min, v_list_min, 'r-', linewidth=2, label='Local Minima', alpha=0.7)
#   plt.plot(t_list_max, v_list_max, 'g-', linewidth=2, label = 'Local Maxima', alpha=0.7)
#   plt.legend()
#   plt.xlabel("Time, t (s)", fontsize=14)
#   plt.ylabel("Displacement, x (m)", fontsize=14)
#   plt.tight_layout(rect=[0.02, 0, 0.98, 1], pad=0.5)
#   plt.show()

  # Plot x,y trajectory
  plt.figure(figsize=(8.5,6))
  y_list = [val[2] for val in V_list]
  plt.plot(x_list, y_list, 'black', linewidth=3, label='Trajectory')
  plt.plot(x_list[0], y_list[0], 'bo', markersize=10, label='Start')
  plt.plot(x_list[-1], y_list[-1], 'ro', markersize=10, label='End')
  plt.xlabel("x displacement (m)", fontsize=18.5)
  plt.ylabel("y displacement (m)", fontsize=18.5)
  plt.tick_params(axis='both', which='major', labelsize=15)
  plt.axis('equal')
  plt.grid(True, alpha=0.3)
  plt.legend(loc=1, prop={'size': 16})
  plt.tight_layout(rect=[0.02, 0, 0.98, 1], pad=0.5)
  plt.show()

  # Print the average value of a period of oscillation
  print("T=",ODE_RK4.period(t_list_min))