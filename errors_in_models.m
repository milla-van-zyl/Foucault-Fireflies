% Define parameters
g = 9.81; % acceleration due to gravity (m/s^2)
L = 67;    % length (m)
m = 28;    % mass (kg)
sig = deg2rad(48.01014); % sigma = co-latitude

% SMALL ANGLE / LINEAR MODEL ==============================
function dRdt = foucaultsmall(R, t, L, g)
% R = [x; y; xdot; ydot];
x = R(1);
y = R(2);
xdot = R(3);
ydot = R(4); 

% Small amplitude
xdotdot = -(g/L)*x;
% Centrifugal force neglected : (omg)^2*R_E*cos(sig)*sin(sig);
ydotdot = -(g/L)*y;

dRdt = [xdot; ydot; xdotdot; ydotdot];
end

% LARGE ANGLE / NONLINEAR MODEL ==============================
function dRdt = foucaultlarge(R, t, L, g)
% R = [x; y; xdot; ydot];
x = R(1);
y = R(2);
xdot = R(3);
ydot = R(4); 

D = (L^2*(xdot^2 + ydot^2) - (x*ydot - y*xdot)^2)/(L^2 - x^2 - y^2);

xdotdot = -(x/L^2)*D - (g/L^2)*x*sqrt(L^2 - x^2 - y^2);
% Centrifugal force neglected : (omg)^2*R_E*cos(sig)*sin(sig);
ydotdot = -(y/L^2)*D - (g/L^2)*y*sqrt(L^2 - x^2 - y^2);

dRdt = [xdot; ydot; xdotdot; ydotdot];
end

% Calculating errors between solutions for displacements 1 to 10m
for i = 1 : 10
    R0 = [i; 0; 0; 0]; % Initial position and velocity
    tspan = linspace(0,5*60,50000); % Time span for simulation (5 mins)
    
    % Solve the coupled equations using ODE solver
    opts = odeset('RelTol', 1e-12,'AbsTol',1e-12);
    [tS, RS] = ode45(@(t, R) foucaultsmall(R, t, L, g), tspan, R0, opts);
    [tL, RL] = ode45(@(t, R) foucaultlarge(R, t, L, g), tspan, R0, opts);
    
    % Extract the results for x and y
    xs = RS(:, 1);
    ys = RS(:, 2);
        
    xL = RL(:, 1);
    yL = RL(:, 2);
    
    err = abs(xs - xL); % Difference between x-solutions
    Emax = max(err) % Maximum difference for each set of solutions
    Emax_percent = Emax*100/max(abs(xL)) % Error as a percentage
end

