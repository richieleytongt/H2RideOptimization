clc, clear

%This script inputs the mass and tire stiffness for a quarter car
%suspension model and uses the interior-point method within fmincon to find
%the optimal stiffnesses of the linear spring and damper connecting the two
%masses.

%The cost function to minimize is J = trace(q*C*P*C') where q is the
%weighing matrix, C is the output matrix, and P is the matrix of unknowns
%in the Lyapunov equation.

%To run this code, MATLAB's Optimization toolbox is necessary.

%This code is calibrated with Example 6.3 in Dynamics and Optimal Control
%of Road Vehicles (Limebeer and Massaro, 2018)

%% Optimization routine
% Fixed parameters
m__s = 120;    %(kg)
m__u = 15;     %(kg)
k__t = 150000; %(N/m)

% Initial guess [k__s, c__s] (N/m) & (Ns/m)
x0 = [15000; 1000];

% Lower bounds (positive stiffness & damping)
lower = [1000; 700];
upper = [];

% weights
q(1) = 70;   %workspace weight
q(2) = 70e5; %roadholding weight

% Optimization options
options = optimoptions('fmincon','Display','iter');

% Run optimizer
[x_opt, J_opt] = fmincon(@(x) H2_cost(x, m__s, m__u, k__t, q), x0, [], [], [], [], lower, upper,  [], options);

k__s_opt = x_opt(1)
c__s_opt = x_opt(2)
J_opt

%% Fermat's theorem

%Taking the partial derivative of the symbolic cost function with respect
%to the spring stiffness yields the location of the stationary point.
%The second derivative being positive allows the below equation to hold for
%the optimal stiffness. The result of the two procedures should be the
%same.

eq652 = (q(2)*k__t*m__s*m__u)/(k__t^2 + q(2)*(m__s + m__u)^2) %Equation 6.52 from the book