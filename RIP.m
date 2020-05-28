%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Initialization
%   First run this before any .slx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all; clc;
%% Loading references and setting parameters
n_states = 4;
n_inputs = 1;
n_outputs = 2;

%% Model Parameters
Rm = 2.6;
Kt = 0.00767;
Km = 0.00767;
Kg = 70;
nu_m = 0.69;
nu_g = 0.9;
Beq = 0.004;
g = 9.81;
Jeq = 0.0044;
r = 0.145;
L = 0.305;
m = 0.210;

a = Jeq + m*r^2;
b = m*L*r;
c = 4*m*L^2/3;
d = m*g*L;
E = a*c - b^2;
G = nu_m*nu_g*Kt*Km*Kg^2/Rm + Beq;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Model and Open Loop Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(n_states);
A(1,3) = 1;
A(2,4) = 1;
A(3,2) = b*d/E;
A(3,3) = -c*G/E;
A(4,2) = a*d/E;
A(4,3) = -b*G/E;

B = zeros(n_states,n_inputs);
B(3) = c*nu_m*nu_g*Kt*Kg/(Rm*E);
B(4) = b*nu_m*nu_g*Kt*Kg/(Rm*E);

C = zeros(n_outputs, n_states);
C(1,1) = 1;
C(2,2) = 1;

D = zeros(n_outputs, n_inputs);
sys = ss(A,B,C,D);

% open loop analysis

pole(sys)
% System has 1 unstable pole, and one marginally stable one.

tzero(sys)
% System has NO transmission zeros

CO = ctrb(A,B);
rank(CO)
% System is controllable 

OB = obsv(A,C);
rank(OB)
% System is observable
% 
% %pzmap + make markers larger
% hpz = pzplot(sys)
% set(hpz.allaxes.Children(1).Children, 'MarkerSize', 12, 'LineWidth', 4)

% % possibly bode diagram ? 
% [MAG,PHASE,W] = bode(sys);
% MM = squeeze(MAG); PP = squeeze(PHASE);
% figure(1); clf(1) 
% loglog(W, MM(1,:), 'b', 'LineWidth', 2);
% hold on;
% loglog(W, MM(2,:), 'g--', 'LineWidth', 2);
% xlabel('Frequency'); ylabel('Magnitude');
% grid on;
% figure(2); clf(2);
% semilogx(W, PP(1,:), 'b', 'LineWidth', 2);
% hold on;
% semilogx(W, PP(2,:), 'g--', 'LineWidth', 2);
% xlabel('Frequency'); ylabel('Phase');
% grid on;
% -> System is said to be minimal as its observable and controllable
% -> System is also detectable and stabilizable as the unstable modes are
%         observable and controllable.



% Control goals : 
%       * we want disturbance rejection of the parameter concerning the rod
%       angle
%
%       * We want to track the setpoint for the parameter conerning the arm
%       angle
%
%       * At the same time we need the rod angle to have a fast reponse
%       such that it is still stabilizatble.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Designing LQR controller
%       entire statevector is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Can we use integral action??

% Constructing the Augmented system
NA = [ 0, C(1,:); zeros(n_states,1), A];

NB = [ D(1,:);B];
sysAG = ss(NA,NB,[0,1,0,0,0;  0,0,1,0,0], zeros(2,1));
CO = ctrb(NA,NB);
rank(CO)
% Rank of the controlability matrix is the same size as the nummber of
% augmented states -> doable

R_i = 1.5;
% %_% Case 1
% Q_i = diag([4 20 0 0]);
% Q_i = [ 1, zeros(1,4); zeros(4,1), Q_i];
% 
% %_% Case 2
% Q_i = diag([4 20 0 0]);
% Q_i = [ 8, zeros(1,4); zeros(4,1), Q_i];

%_% Case 3
Q_i = diag([4 20 0 600]);
Q_i = [ 8, zeros(1,4); zeros(4,1), Q_i]; 

[K_d_Int ,~, CLP] = lqr(NA,NB,Q_i, R_i);

Ki = K_d_Int(:,1);
Ks = K_d_Int(:,2:end);

%close loop poles
K_adj = -NB*K_d_Int; 
Closed_A = NA + K_adj;
polesOfClosedLoopSystem = eig(Closed_A) % --> poles of closed loop system

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Designing LQR controller
%       State vector should be calculated using backward differences
%       no estimation needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc =2;% 0.4;%10;% 0.2; %1.5
wc = 2*pi*fc;
Ts = 1/200;