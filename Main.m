clc 
clear 
close all

%% Simulation setup
TimeSpan = 10;
DT = 1e-3;
Plant = 1;

%% Declaration of the vector dimensions
n = 5; %states
p  = 1; %input u
q  = 3; %output y
m = 1;  %error e 
nd = 2; %disturbs d
r = nd + q + m; %exogenous w


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
A_lift = 3; %m
Cd_lift = 0.5; 
rho_lift = 1.225; %kg/m^3
%speed of the car
v = 10; %[m/s]


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
nu0 = [0;0;0];
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
plant0 = []
%% LINEARIZATION POINT
x_tilde_init = x_init - x0;
u_tilde = 0;
y_tilde = 0; 

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
%% MATRICES OF e - e(1) = x(1) + vi(1) - r_l

Ce = [1 0 0 0 0];  

De1 = [0];

De2 = [0 0 1 0 0 -1];

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

%% STATE FEEDBACK CONTROL + INTEGRAL ACTION
Ae = [A zeros(n,m);
      Ce zeros(m,m)];  % dot(x,eta) = Ae*(x,eta)

Be = [B1
       zeros(m,p)]; %dot(x,eta) = Be*u

% eps=C_eps*[x,eta]
% eps=[x1,x2,x3,x4,x5,eta1,Ddot(z_s),z_s-z_r=x1+x3]
Ceps = [1 0 0 0 0 0
        0 1 0 0 0 0
        0 0 1 0 0 0
        0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        -ks/ms -betas/ms 0 0 Ap/mi 0
        1 0 1 0 0 0];

% eps=D1eps*u
D1eps = [0; 0; 0; 0; 0; 0; 0; 0];

% DRIVE MODES
drive_mode = 0; % 0 = comfort; 1 = off-road; 2 = race
switch drive_mode 
    % not to penalize eps, we put it equal to 1e4
    case 0
        eps1max = 1e4; % susp. deflection
        eps2max = 1; % susp. speed
        eps3max = 1e4; % tire deflection
        eps4max = 0.1; % tire deflection speed PEN
        eps5max = 1e6; % actuator force (do not pen)
        eps6max = 0.1; % integral of the position error PEN 
        eps7max = 0.01*g; % sprung mass acceleration PEN
        eps8max = 1; % sprung mass height 
     
    case 1
        eps1max = 1e4; % susp. deflection
        eps2max = 1e4; % susp. speed
        eps3max = 0.001; % tire deflection PEN 
        eps4max = 0.001; % tire speed PEN 
        eps5max = 1e6; % actuator force (do not pene)
        eps6max = 0.001; % integral of the position error PEN
        eps7max = 1e4; % sprung mass acceleration
        eps8max = 1e4; % sprung mass height 

     case 2
        eps1max = 0.001; % susp. deflection PEN
        eps2max = 0.001; % susp. speed PEN
        eps3max = 1e5; % tire deflection 
        eps4max = 1e5; % tire speed
        eps5max = 1e6; % actuator force (do not pene)
        eps6max = 0.001; % integral of the position error PEN
        eps7max = 1e4; % sprung mass acceleration
        eps8max = 0.001; % sprung mass height PEN
end

Q = inv(8*diag([eps1max^2,eps2max^2,eps3max^2,eps4max^2,eps5max^2, ...
    eps6max^2,eps7max^2,eps8max^2]));

umax = 0.01; % source: Road Adaptive....
R = inv(umax^2);
barR = R+D1eps.'*Q*D1eps;

alpha = 0; %we already have RE(lambda)>0

Am = Ae+alpha*eye(n+m);
Em = eye(n+m);
Bm = Be;
Gm = 0;
Qm = Ceps.'*Q*Ceps;
Sm = Ceps.'*Q*D1eps;
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

% DA FINIRE con i dati dei sensori

w1max = 500; %maximum road acceleration
w2max = 70; %maximum lift force (fixed)

std_pot = 0.001; % [m] potentiomenter standard deviation
std_laser = 0.001; %standard deviation laser
std_acc = 110e-6*sqrt(200); % [m/s^2] accelerometer standard deviation

Qd = diag([w1max^2,w2max^2,0,0,0,0]);
Rd = diag([std_pot^2,std_laser^2,std_acc^2]);
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