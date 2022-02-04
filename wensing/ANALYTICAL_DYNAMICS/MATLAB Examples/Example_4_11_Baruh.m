%% Example 4.11 - Baruh
% Code Suplement to find equations of motion using Lagrange's Equations.
clc;
clear all;

%% Setup:
% Symbolic parameters:
syms M m k c b g real
syms x   theta   real
syms dx  dtheta  real
syms ddx ddtheta real
syms F Fx real

% Symbolic groupings:

q   = [x theta]';       % Generalized coordinates
dq  = [dx dtheta]';     % Velocity
ddq = [ddx ddtheta]';   % Acceleration
F_ext = [F 0 0]';       % External force
F_d  = [Fx 0 0]';       % Force from dashpot

p = [M; m; k; c; b; g]; % System Parameters

%% Handy anonymous functions:
ddt = @(x)( jacobian(x,q)*dq + jacobian(x,dq)*ddq); % Function to compute time derivatives
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));          % Force contributions to generalized forces
                                                    % Equation [4.6.3] - Baruh
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));         % Moment contributions to generalized forces

%% Kinematics:
% ** Variables expressed with respect to an inertial frame O: X/Y/Z (Cartesian Coordiantes).  
O_ehat1 = [sin(theta) -cos(theta) 0]'; % Unit vector - Position of the CoM of the bar
O_ghat = [0 -1 0]';                    % Unit vector - Gravity

O_r_cart = [x 0 0]';                   % Position of the cart
O_v_cart = ddt(O_r_cart);                % Velocity of the cart

O_r_bar = O_r_cart + O_ehat1*b/2;          % Position CoM of the bar
O_v_bar = ddt(O_r_bar);                  % Velocity CoM of the bar

O_r_tip = O_r_cart + O_ehat1*b;            % Position tip of the bar

%% Lagrangian:
% Kinetic and Potential Energy of the cart:
T_cart = simplify(1/2 * M * (O_v_cart.' * O_v_cart));
V_cart = simplify(1/2 * k * (O_r_cart.' * O_r_cart));

% Kinetic and Potential Energy of the bar:
T_bar = simplify(1/2 * (1/12 * m * b^2) * dtheta^2 + 1/2 * m * (O_v_bar.' * O_v_bar)); 
V_bar = simplify(m * g * dot(O_r_bar,-O_ghat));

% Compute Lagrangian and Total Energy:
L = simplify(T_cart + T_bar - V_cart - V_bar); % Lagrangian
E = simplify(T_cart + T_bar + V_cart + V_bar); % Energy

%% Generalized Forces:
O_Fxhat = [-1 0 0]';                        % Unit vector - Dashpot Force
O_Fx2 = simplify(c * dot(O_v_bar,O_Fxhat)); % Force due to dashpot
O_Fd2 = subs(F_d,Fx,O_Fx2);                 % Force substituted

Q_F = F2Q(F_ext,O_r_tip);     % Generalized forces due to external force
Q_Fd = F2Q(O_Fd2,O_r_bar);    % Generalized forces due to dashpot force

% ** There is a mistake in Baruh with respect to the generalized forces.
% ** You need to check that the force units agree.
% ** Take into account that the dashpot force is given in terms of the
% ** velocity of the CoM of the bar, which has an angular-change component.

Q = Q_F + Q_Fd;             % Total generalized forces

%% Dynamics:
dL_dq = jacobian(L,q)';  % Lagrangian component 1
dL_dqd= jacobian(L,dq)'; % Lagrangian component 2

%% Equations of motion:
eom = ddt(dL_dqd) - dL_dq - Q






