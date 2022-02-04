% Symbolic math toolbox example to support Lecture 2

clear
syms M m L theta theta_dot theta_ddot x x_dot x_ddot F T N g real

q = [x ; theta]; % Degrees of freedom. What we shall later call generalized coordinates
q_dot = [x_dot ; theta_dot];
q_ddot = [x_ddot ; theta_ddot];

% Standard Unit Vectors
i_hat = [1 ; 0];
j_hat = [0 ; 1];

% Handy unit vectors 
e_T    = cos(theta)*j_hat - sin(theta)*i_hat; % Unit vector from mass to cart
e_perp = cos(theta)*i_hat + sin(theta)*j_hat; % Perpendicular to it

% Autonomous function for taking time derivatives.
ddt = @(f) jacobian(f,q)*q_dot + jacobian(f,q_dot)*q_ddot;

p_cart = [x ; 0];
p_mass = p_cart - L* e_T;

v_cart = ddt(p_cart);
v_mass = ddt(p_mass);

a_cart = ddt(v_cart);
a_mass = ddt(v_mass);

F_net_cart = N*j_hat - M*g*j_hat + F*i_hat - T*e_T;
F_net_mass = T*e_T - m*g*j_hat;

% First equation of motion from class, add the x-direction equations of motion for cart and mass
EOM1 = simplify( dot(i_hat, M*a_cart + m*a_mass)) == simplify( dot(i_hat, F_net_mass+F_net_cart))

% Second equation of motion from class, consider Newton's law for the second mass in the e_perp direction
EOM2 = simplify( dot(e_perp, m*a_mass) ) == simplify( dot(e_perp, F_net_mass) )
