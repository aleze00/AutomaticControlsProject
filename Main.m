clc 
clear all
close all

%% Simulation setup
TimeSpan = 10;
DT = 1e-3;
Plant = 1;

%% Declaration of the vector dimensions
n = 4; 
p  = 1;
q  = 2;
m = 3;
nd = 2;
r = nd + q + m;


%% Initial conditions
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
%% Plant Parameters
ks = 29300; %spring coefficient
kt = 290000; %tyre elastic coefficient
mu = 38;
ms = 395;
l0s = 0.5;
l0t = 0.10;
g = 9.81;
betas = 3000; %damping coefficient 

%% Linearization initial conditions
% in order to compute the equilibrium, we put u=0, d=0, dot{x} = 0 in the
% non-linearized equation, find x0.
u0 = 0;
d0 = [0;0];
nu0 = [0;0];
r0 = l0s-g*ms/ks;
w0 = [d0; nu0; r0];

%obtained from paperwork
x0 = [l0s-g*ms/ks
    0
    l0t-g*(ms+mu)/kt
    0];
%from the model equations
y0 = [x0(1)
    x0(1)+x0(3)
    (1/ms)*(d0(2)-ms*g-ks*(x0(1)-l0s)-betas*x0(2)+u0(1));
    ];


