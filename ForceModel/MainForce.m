clc 
clear 
close all

%% Simulation setup
TimeSpan = 10;
DT = 1e-4;
Plant = 1;


%% Declaration of the vector dimensions
n = 5;                                      % DA RIFARE CON I NUOVI STATI
p  = 1;                                     % input u
q  = 3;                                     % output y
m = 1;                                      % error e 
nd = 2;                                     % disturbs d
r = nd + q + m;                             % exogenous w


%% Initial conditions
zu_init = 0.342;
zs_init = 0.641 + zu_init;
vs_init = 0;
vu_init = 0;
zr_init = 0;
dotzr_init = 0;
Force_init = 0;

x_init = [zs_init-zu_init
    vs_init-vu_init
    zu_init-zr_init
    vu_init-dotzr_init];

%% Plant Parameters
ks = 50000; % spring coefficient
kt = 500000; % tyre elastic coefficient
mu = 100;
ms = 300;
l0s = 0.7; % [m] lenght of the sprung for which the force is 0
l0t = 0.35; % [m]
g = 9.81;
betas = 500; % damping coefficient 

A_lift = 3;
Cd_lift = 0.5; 
rho_lift = 1.225; % kg/m^3
v = 10; % [m/s] speed of the car
Lift = 0.5*Cd_lift*(v^2)*rho_lift;

%% Linearization initial conditions

% in order to compute the equilibrium, we put u=0, d=0, dot{x} = 0 in the
% non-linearized equation, find x0.

u0 = 0;
d0 = [0;Lift];
nu0 = [0;0;0];
r0 = l0s-g*ms/ks; % lenght of the sprung when the car load is stationary (same as the equilibrium point x0(1)
w0 = [d0; nu0; r0];

%obtained from paperwork DA RIFARE CON I NUOVI STATI
x0 = [l0s-g*ms/ks
    0
    l0t-g*(ms+mu)/kt
    0];

%from the model equations
y0 = [x0(1)                                 % suspension lenght (potentiometer)
    x0(1)+x0(3)                             % car height (laser)
    (1/ms)*(Lift-ms*g-ks*(x0(1)-l0s))];     % passenger acceleration

e0 = y0(1) - r0;  
   
    
%% LINEARIZATION POINT
x_tilde_init = x_init - x0;
u_tilde = 0;
y_tilde = 0; 

%% LINEARISED PLANT

barm = (mu+ms)/ms;

% DA RIFARE CON I NUOVI STATI

A = [0 1 0 0
    -(ks*(ms+mu))/(ms*mu) -(betas*(ms+mu))/(ms*mu) kt/mu 0
    0 0 0 1
    ks/mu betas/mu -kt/mu 0]; 

B1 = [0;0;0;0];

B2 = [0 0 0 0 0 0
    0 1/ms 0 0 0 (ks*(ms+mu))/(ms*mu)
    0 0 0 0 0 0
    -1 0 0 0 0 -ks/mu];

C = [1 0 0 0 0
    1 0 1 0 0
    -ks/ms -betas/ms 0 0 1/ms];

D1 = [0;0;0];

D2 = [0 0 1 0 0 0
    0 0 0 1 0 0
    0 1/ms 0 0 1 ks/ms];

%% MATRICES OF e - e(1) = x(1) + vi(1) - r_l
Ce = [1 0 0 0 0];  

De1 = [0];

De2 = [0 0 1 0 0 -1];

%% REACHABILITY CHECK for A,B1
 R = ctrb(A,B1);
 if rank(R)==length(R)
     disp("FULLY REACHABLE")
 else
    rankR = rank(R);
    disp('NOT FULLY REACHABLE')
    imR = orth(R);              %pg.29 PDF toni
    imRorth = null(R.');         %ker(R)
    invTR = [imR imRorth];      %T_R^-1
    barA = inv(invTR)*A*invTR;  
    A22 = barA(rankR+1:end, rankR+1:end);  %non-reachable part of barA
    eA22 = eig(A22);
    for i = 1:length(eA22)  
        if real(eA22(i)) >= 0
            disp('NOT STABILISABLE')
        end
    end
 end

%% OBSERVABILTY CHECK for A,C
O = obsv(A,C);
rankO = rank(O);
if rankO == length(A)
    disp('FULLY OBSERVABLE')
else
    disp('NOT FULLY OBSERVABLE')
end

%% STATE FEEDBACK CONTROL + INTEGRAL ACTION
Ae = [A zeros(n,m);
      Ce zeros(m,m)];

Be = [B1
       zeros(m,p)];

Ceps = [1 0 0 0 0 0
        0 1 0 0 0 0
        0 0 1 0 0 0
        0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        -ks/ms -betas/ms 0 0 1 0
        1 0 1 0 0 0];

D1eps = [0; 0; 0; 0; 0; 0; 0; 0];

Q = inv(8*diag([eps1max^2,eps2max^2,eps3max^2,eps4max^2,eps5max^2, ...
    eps6max^2,eps7max^2,eps8max^2]));

umax = 0.01; % source: Road Adaptive....
R = inv(umax^2);
barR = R+D1eps.'*Q*D1eps;
alpha = 0; 
Am = Ae+alpha*eye(n+m);
Em = eye(n+m);
Bm = Be;
Gm = 0;
Qm = Ceps.'*Q*Ceps;
Sm = Ceps.'*Q*D1eps;
Rm = barR;

%% REACHABILITY CHECK of Am,Bm
 lengthAm = length(Am);
 R = ctrb(Am,Bm);
 rankR = rank(R);
 disp('---REACHABILITY CHECK of Am,Bm---');
 if rank(R)==lengthAm
     disp("FULLY REACHABLE")
 else
    disp('NOT FULLY REACHABLE')
    imR = orth(R);              
    imRorth = null(R.');         %ker(R)
    invTR = [imR imRorth];      %T_R^-1
    barA = inv(invTR)*Am*invTR;  
    A22 = barA(rankR+1:end, rankR+1:end);  %non-reachable part of barA
    eA22 = eig(A22);
    for i = 1:length(eA22)  
        if real(eA22(i)) >= 0
            disp('-NOT STABILISABLE thanks to eigenvalue number: ');
            disp(i);
        else disp('STABILISABLE');
        end
    end
 end

%% OBSERVABILTY CHECK of Am, Ceps
O = obsv(Am,Ceps);
rankO = rank(O);
kerO = null(O);
disp('---OBSERVABILTY CHECK of Am, Ceps---');
if rankO == length(Am)
    disp('FULLY OBSERVABLE')
else
    disp('NOT FULLY OBSERVABLE')
    A11 = barA(1:(length(Ceps)-rankO), 1:(length(Ceps)-rankO));  %non-observable part of barA
    eA11= eig(A11);
    for i = 1:length(eA11)  
        if real(eA11(i)) >= 0
            disp('-NOT DETECTABLE thanks to eigenvalue number: ');
            disp(i);
        end
    end
end

%% SOLVING RICCATI EQUATION
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