function J = H2_cost(x, m__s, m__u, k__t, weight)

k__s = x(1);
c__s = x(2);

%Input C Matrix
C = [-k__s/m__s, -c__s/m__s, 0, c__s/m__s; 1, 0, 0, 0; 0, 0, 1, 0];

%Input the elements of the P matrix. [solved symbolically in Maple]
m = m__s + m__u; %computational efficiency
P11 = m/(2*c__s);
P12 = -m/(2*m__s);
P13 = (k__s*m^2 - m__s*m__u*k__t)/(2*m__s*k__t*c__s);
P14 = -m/(2*m__s);
P22 = (k__s^2*m^2 + k__t*(m*c__s^2 + m__s^2*k__s))/(2*m__s^2*k__t*c__s);
P23 = m__u/(2*m__s);
P24 = (k__s^2*m^2 + k__t*(m*c__s^2 - m__s*m__u*k__s))/(2*m__s^2*k__t*c__s);
P33 = (k__s^2*m^3 + k__t*c__s^2*m^2 + m__s*m__u*k__t*(m__s*k__t  - 2*m*k__s))/(2*m__s^2*k__t^2*c__s);
P34 = -0.5;
P44 = (k__s^2*m^2 + k__t*(c__s^2*m + m__s^2*(k__t - k__s) - 2*m__s*m__u*k__s))/(2*m__s^2*k__t*c__s);

% Construct the P matrix
P = [P11, P12, P13, P14; P12, P22, P23, P24; P13, P23, P33, P34; P14, P24, P34, P44];

% Construct the q matrix
q = diag([1 weight(1) weight(2)]);

%Cost function to minimize
J = trace(q*C*P*C.');