% Define parameters
g = 9.81; % acceleration due to gravity (m/s^2)
L = 67;    % length (m)
% gam = damping coefficient
m = 28;    % mass (kg)
omg = 2*pi/86400; % angular frequency of Earth's rotation
sig = deg2rad(48.01014); % sigma = co-latitude

function dRdt = foucault(R, t, omg, sig, L, g)
% R = [x; y; xdot; ydot];
x = R(1);
y = R(2);
xdot = R(3);
ydot = R(4); 

% Small amplitude
xdotdot = 2*omg*cos(sig)*ydot - (x/L)*g;
% Centrifugal force neglected : (omg)^2*R_E*cos(sig)*sin(sig);
ydotdot = -2*omg*cos(sig)*xdot - (y/L)*g;

dRdt = [xdot; ydot; xdotdot; ydotdot];
end

% Initial conditions
R0 = [1; 0; 0; 0]; % Initial position and velocity
T = 2*pi* sqrt(L/g); % Time period
tspan = 0 : T : 200000; % Time span for simulation, plotting one point per period

% Solve the coupled equations using ODE solver
opts = odeset('RelTol', 1e-12,'AbsTol',1e-12);
[t, R] = ode45(@(t, R) foucault(R, t, omg, sig, L, g), tspan, R0, opts);

% Extract the results for x and y
x = R(:, 1);
y = R(:, 2);

% Plot the results
figure;
plot(tspan/3600, x, '-', 'LineWidth', 1.5);
hold on;
plot(tspan/3600, y, '-', 'LineWidth', 1.5);
hold on;
Y = sin(omg*tspan);
plot(tspan/3600, Y, '-', 'LineWidth', 1.5);
xlabel('Time in hours', 'Interpreter','latex');
ylabel('Displacement', 'Interpreter','latex');
title('Foucault Pendulum Motion', 'Interpreter','latex');
legend('Pendulum precession, x', 'Pendulum precession, y', 'Earth rotation','Interpreter','latex','Location', 'Best');
fontsize(16,"points")
xlim([0 55]);
grid on;