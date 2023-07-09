clc
close all
clear all

%% SIMULATION SETUP

TimeSpan = 10; % [s] simulation time span
DT = 1e-3; % [s] sample time fixed-step integration
PLANT = 1; % 0 = linear, 1 = nonlinear

%% DECLARATION OF VECTOR DIMENSIONS

n = 4; % state
p = 1; % control
q = 2; % measurement
m = 1; % regulated output
nd = 1; % disturbance
r = nd+q+m; % exogenous

%% INITIAL CONDITIONS FOR THE NONLINEAR PLANT

zu_init = 0.10;
zs_init = 0.5+zu_init;
vs_init = 0;
vu_init = 0;
zr_init = 0;
dotzr_init = 0;

x_init = [zs_init-zu_init
    vs_init-vu_init
    zu_init-zr_init
    vu_init-dotzr_init];

% x_init = [0.3467
%     0
%     0.0819
%     0];

%% NOMINAL PLANT PARAMETERS USED TO DESIGN THE CONTROL SYSTEM

ks = 29300;
kt = 290000;
mu = 38;
ms = 395;
ells = 0.5;
ellt = 0.10;
g = 9.81;
bs = 3000;


%% LINEARISATION CONDITIONS

u0 = 0;
d0 = 0;
nu0 = [0;0];
r0 = ells-g*ms/ks;
w0 = [d0; nu0; r0];
x0 = [ells-g*ms/ks
    0
    ellt-g*(ms+mu)/kt
    0];
y0 = [x0(1)
    -g-ks/ms*(x0(1)-ells)-bs/ms*x0(2)+u0/ms+g];
e0 = x0(1)-r0;


%% LINEARISED PLANT

% state
a21 = -(mu+ms)/(ms*mu)*ks;
a22 = -(mu+ms)/(ms*mu)*bs;
a23 = kt/mu;
a41 = ks/mu;
a42 = bs/mu;
a43 = -kt/mu;

A = [0 1 0 0;
    a21 a22 a23 0;
    0 0 0 1;
    a41 a42 a43 0];
[V,Vn,J] = JCF(A)

b12 = (mu+ms)/(ms*mu);
b14 = -1/mu;
b24 = -1;

B1 = [0; b12; 0; b14];

B2 = [0 0 0 0;
    0 0 0 0;
    0 0 0 0;
    b24 0 0 0];

% output
C = [1 0 0 0;
    -ks/ms -bs/ms 0 0];
D1 = [0; 1/ms];
D2 = [0 1 0 0
    0 0 1 0];

% error
Ce = [ 1 0 0 0];
De1 = 0;
De2 = [0 0 1 -1];

% initial conditions
tilde_x_init = x_init-x0;

%% STATE FEEDBACK CONTROL + INTEGRAL ACTION

Ae = [A zeros(n,m)
    Ce zeros(m,m)];

Be = [B1;
    zeros(m,p)];

Ceps = [-ks/ms -bs/ms 0 0 0;
           1 0 0 0 0;
           0 0 1 0 0;
           0 0 0 0 1;
           0 0 0 1 0;
           0 1 0 0 0];
%    1 0 1 0 0;
%    0 1 0 1 0];

Deps =  [1/ms; 0; 0; 0;0; 0]; %0;0];

scenario = 0; % 0 = city; 1 = race; 2 = off roads
switch scenario
    case 0
        eps1max = 0.02*g; % [m/s^2] sprung mass acceleration max all.
        eps2max = 1e4; % [m] suspension deflection max all.
        eps3max = 1e4; % [m] tire deflection max all.
        eps4max = 0.1; % [m*s] integral of the position error max all.
        eps5max = 0.01; % [m/s] tire variation speed max all.
        eps6max = 1; % to penalise the susp. speed
%        eps7max = 0.01;
    case 1
        eps1max = 0.1*g; % [m/s^2] sprung mass acceleration max all.
        eps2max = 1e-4; % [m] suspension deflection max all.
        eps3max = 1e5; % [m] tire deflection max all.
        eps4max = 1e-4; % [m*s] integral of the position error max all.
        eps5max = 0.005; % [m/s] tire variation speed max all.
%         eps6max = 1e4;
%         eps7max = 1e4;
    case 2
        eps1max = 1e4; % [m/s^2] sprung mass acceleration max all.
        eps2max = 1e4; % [m] suspension deflection max all.
        eps3max = 0.001; % [m] tire deflection max all.
        eps4max = 0.001; % [m*s] integral of the position error max all.
        eps5max = 0.001; % [m/s] tire variation speed max all.
%         eps6max = 1e4;
%         eps7max = 1e4;
end

Q = inv(6*diag([eps1max^2,eps2max^2,eps3max^2,eps4max^2,eps5max^2,eps6max^2]));%,eps6max^2,eps7max^2]));

umax = 10000; % [N] active Force max all.
R = inv(umax^2);

barR = R+Deps.'*Q*Deps;

alpha = 0;

Am = Ae+alpha*eye(n+m);
Em = eye(n+m);
Bm = Be;
Gm = 0;
Qm = Ceps.'*Q*Ceps;
Sm = Ceps.'*Q*Deps;
Rm = barR;

[X,Km,L] = icare(Am,Bm,Qm,Rm,Sm,Em,Gm);
K = -Km;
KS = K(:,1:n);
KI = K(:,n+1:n+m);

%% OBSERVER

lambda_d = 0;
Ad = A.';
Bd = C.';
Cd = B2.';
Dd = D2.';

w1max = 100;
std_pot = 0.001; % [m] potentiomenter standard deviation
std_acc = 0.05*g; % [m/s^2] accelerometer standard deviation

Qd = diag([w1max^2,0,0,0]);
Rd = diag([std_pot^2,std_acc^2]);
barRd = Rd + Dd.'*Qd*Dd;

Am = Ad+lambda_d*eye(n);
Em = eye(n);
Bm = Bd;
Gm = 0;
Qm = Cd.'*Qd*Cd;
Sm = Cd.'*Qd*Dd;
Rm = barRd;
[X,Km,L] = icare(Am,Bm,Qm,Rm,Sm,Em,Gm);
KO = Km.';

AO = A-KO*C;
BO = B1-KO*D1;
CO = eye(n);
DO = zeros(n,q);

% initial condition
% xO_init = zeros(n,1);
xO_init = tilde_x_init;