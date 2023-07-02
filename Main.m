clc 
clear 
close all

%% Simulation setup
TimeSpan = 10;
DT = 1e-3;
Plant = 1;

%% Declaration of the vector dimensions
n = 5; 
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
Pl_init = 0;

x_init = [zs_init-zu_init
    vs_init-vu_init
    zu_init-zr_init
    vu_init-dotzr_init
    Pl_init];

%% Plant Parameters
ks = 29300; %spring coefficient
kt = 290000; %tyre elastic coefficient
mu = 38;
ms = 395;
l0s = 0.5;
l0t = 0.10;
g = 9.81;
betas = 3000; %damping coefficient 
mi = 1e-7; % scale coefficient to improve numerical conditioning of P
alfa = 4.515e13; % N/m^5 

%%ACTUATOR PARAMETERS
 beta = 1; %1/sec
 gamma = 1.545e9; % N/m^(5/2)kg^(1/2)
 Ap = 3.35e-4; % m^2
 Ps = 10342500; %Pa
 rho = 865; % Kg/m^3

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
    0
    0];
%from the model equations
y0 = [x0(1)
    x0(1)+x0(3)
    (1/ms)*(d0(2)-ms*g-ks*(x0(1)-l0s)-betas*x0(2)+u0(1));
    ];

e0 = [y0(3)
    y0(1) - r0      
    y0(2)-y0(1)
    ];

%% LINEARISED PLANT

A = [0 1 0 0 0
    -(ks*(ms+mu))/(ms*mu) -(betas*(ms+mu))/(ms*mu) kt/mu 0 -(Ap/mi*(ms+mu))/(ms*mu)
    0 0 0 1 0
    ks/mu betas/mu -kt/mu 0 -Ap/(mu*mi)
    0 -mi*alfa*Ap 0 mi*alfa*Ap -beta];

B1 = [0;0;0;0;mi*gamma*sqrt(Ps/rho)];

B2 = [0 0 0 0 0 0
    0 1/ms 0 0 0 (ks*(ms+mu))/(ms*mu)
    0 0 0 0 0 0
    -1 0 0 0 0 -ks/mu
    0 0 0 0 0 0];

C = [1 0 0 0 0
    1 0 1 0 0
    -ks/ms -betas/ms 0 0 Ap/(mi*ms)];

D1 = [0;0;0];

D2 = [0 0 1 0 0 0
    0 0 0 1 0 0
    0 1/ms 0 0 1 ks/ms];

C_e = [1 0 0 0 0];

D_e1 = [0];

D_e2 = [0 0 1 0 0 -1];
 
%[V,Vn,J] = JCF(A);

%% REACHABILITY CHECK

 R = ctrb(A,B1);
 if rank(R)==length(R)
     disp("FULLY REACHABLE")
 else
    disp('NOT FULLY REACHABLE')
    imR = orth(R);              %pg.29 PDF toni
    imRorth = ker(R.');         %ker(R)
    invTR = [imR imRorth];      %T_R^-1
    barA = inv(invTR)*A*invTR;  
    A22 = barA(rankR+1:end, rankR+1:end);  %non-reachable part of barA
    eA22 = eig(A22);
    %Stabilizabilty Check
    for i = 1:length(eA22)  
        if real(eA22(i)) >= 0
            disp('NOT STABILISABLE')
        end
    end
 end

%% OBSERVABILTY CHECK

O = obsv(A,C);
rankO = rank(O);
if rankO == length(A)
    disp('FULLY OBSERVABLE')
else
    disp('NOT FULLY OBSERVABLE')
end


