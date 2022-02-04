%% Example 4.12 - Baruh
% Code Suplement to find equations of motion using Lagrange's Equations.
clc;
clear all;

%% Setup:
% Symbolic parameters:
syms M m L mu psi g real
syms r   theta   real
syms dr  dtheta  real
syms ddr ddtheta real
syms F Ff real

% Symbolic groupings:
q   = [r theta]';       % Generalized coordinates - Polar coordiantes
dq  = [dr dtheta]';     % Velocity
ddq = [ddr ddtheta]';   % Acceleration

p = [M; m; L; mu; psi; g]; % System Parameters

%% Handy anonymous functions:
ddt = @(x)(jacobian(x,q)*dq + jacobian(x,dq)*ddq);  % Function to compute time derivatives
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));          % Force contributions to generalized forces
                                                    % Equation [4.6.3] - Baruh
F2Qb = @(F,v) simplify(jacobian(v,dq)'*(F));        % Force contributions to generalized forces
                                                    % F2Qb is used in case of a moving frame
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));         % Moment contributions to generalized forces

%% Kinematics:
% ** Variables expressed with respect to moving frame A: r/theta/Z (Polar coordinates).
A_r_collar = [r 0 0]';                                        % Position of the collar
A_w_collar = [0 0 dtheta]';                                   % Angular velocity 
A_v_collar = ddt(A_r_collar) + cross(A_w_collar,A_r_collar);  % Velocity 

A_r_P = [L 0 0]';                                             % Position rod's tip
A_v_P = ddt(A_r_P) + cross(A_w_collar,A_r_P);                 % Velocity 

A_r_rod = [L/2 0 0]';                                         % Position of the rod's CoM

A_ghat = [cos(theta) -sin(theta) 0]';                         % Unit vector - Gravity                                % Unit vector - Gravity
%% Lagrangian:
% Kinetic and Potential Energy of the collar:
T_collar = simplify(1/2 * m * (A_v_collar.' * A_v_collar));
V_collar = simplify(m * g * dot(A_r_collar,-A_ghat));

% Kinetic and Potential Energy of the rod:
T_rod = simplify(1/2 * (1/3 * M * L^2) * dtheta^2); % With respect to O
V_rod = simplify(M * g * dot(A_r_rod,-A_ghat));

% Compute Lagrangian and Total Energy:
L = simplify(T_collar + T_rod - V_collar - V_rod); % Lagrangian
E = simplify(T_collar + T_rod + V_collar + V_rod); % Energy

%% Generalized Forces:
A_F_ext = F*[cos(psi) sin(psi) 0]';     % External Force
A_F_fri  = Ff*[-sign(dr) 0 0]';         % Friction Force

% ** Notice the use of F2Qb below due to the variables expressed in the 
% ** moving frame: 
Q_F = F2Qb(A_F_ext,A_v_P);        % Generalized forces due to external force
Q_Fd = F2Qb(A_F_fri,A_v_collar);  % Generalized forces due to friction force

Q = Q_F + Q_Fd;                   % Total generalized forces

%% Dynamics:
dL_dq = jacobian(L,q)';  % Lagrangian component 1
dL_dqd= jacobian(L,dq)'; % Lagrangian component 2

%% Equations of motion:
eom = ddt(dL_dqd) - dL_dq - Q






