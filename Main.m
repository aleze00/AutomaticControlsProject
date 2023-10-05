clc 
clear 
close all

% più faccio crescere u e più mi allontano dall'equilibrio
% per questo con il gain succede un casino
% forse potremmo alzare gamma?
% quello non lineare dovrebbe migliorare con il gain

%% Simulation setup
TimeSpan = 10;
DT = 1e-4;
Plant = 1;

%% Declaration of the vector dimensions
n = 5; %states
p  = 1; %input u
q  = 3; %output y
m = 1;  %error e 
nd = 2; %disturbs d
r = nd + q + m; %exogenous w

% x = sym('x', [5 1], 'real');
% w = sym('w', [6 1], 'real');
% syms u real

% syms ks kt mu ms l0s l0t g betas Lift real
% syms beta mi gamma Ap Ps alfa rho real

%% Plant Parameters
ks = 29000; % [N/m] spring coefficient
kt = 293000; % [N/m] tyre elastic coefficient
mu = 38;
ms = 295;
l0s = 0.7; % [m] lenght of the sprung for which the force is 0
l0t = 0.35; % [m]
g = 9.81;
betas = 3000; % damping coefficient 

A_lift = 3;
Cd_lift = 0.5; 
rho_lift = 1.225; % kg/m^3
v = 10; % [m/s] speed of the car
Lift = 0.5*A_lift*Cd_lift*(v^2)*rho_lift;

%% ACTUATOR PARAMETERS
beta = 1; %1/sec
mi = 1e-7; % scale coefficient to improve numerical conditioning of P
gamma = 1.545e9; % N/m^(5/2)kg^(1/2)
Ap = 3.35e-4; % m^2
alfa = 4.515e9; % N/m^5

% Ps = 10342500; % Pa
% rho = 865; % Kg/m^3 hydraulic fluid density found on internet
% P = pressure drop across pistons
% x5 = mi*P

% %% State Space Model
% 
% dotxstar1 = x(2);
% dotxstar2 = ((mu+ms)/(ms*mu))*(-ks*(x(1)-l0s)-betas*x(2))-1/mu*(-kt*(x(3)-l0t))+ w(2)/ms + ((mu+ms)/(ms*mu))*(Ap/mi)*x(5);
% dotxstar3 = x(4);
% dotxstar4 = -g + ks/mu*(x(1)-l0s) + betas/mu*x(2) + 1/mu*(-kt*(x(3)-l0t)) - w(1) - 1/mu*(Ap/mi)*x(5);
% dotxstar5 = -beta*x(5) - mi*alfa*Ap*x(2) + mi*gamma*u;
% 
% dot_xstar = [dotxstar1;dotxstar2;dotxstar3;dotxstar4;dotxstar5];
% 
% %% Matrixes Calculation using Numbers
% 
% Astar = jacobian(dot_xstar, x);
% B1star = jacobian(dot_xstar, u);
% B2star = jacobian(dot_xstar, w);

%% Linearization initial conditions
% in order to compute the equilibrium, we put u=0, dot{x} = 0 in the
% non-linearized equation, find x0.
u0 = 0;
d0 = [0;Lift];
nu0 = [0;0;0];
r0 = - g*ms/ks + l0s + Lift/ks; % lenght of the sprung when the car load is stationary (same as the equilibrium point x0(1))
w0 = [d0; nu0; r0];

%obtained from paperwork
x0 = [- g*ms/ks + l0s + Lift/ks
    0
    Lift/kt + l0t - g*(mu+ms)/kt
    0
    0];

%from the model equations
y0 = [x0(1)                                         % suspension lenght (potentiometer)
    x0(1)+x0(3)                                     % car height (laser)
    (1/ms)*(Lift-ms*g-ks*(x0(1)-l0s))];             % passenger acceleration
    
e0 = y0(1) - r0;

%% LINEARISED PLANT

N = (ms+mu)/(ms*mu);

A = [0 1 0 0 0
    -ks*N -betas*N kt/mu 0 Ap/mi*N
    0 0 0 1 0
    ks/mu betas/mu -kt/mu 0 -Ap/(mu*mi)
    0 -mi*alfa*Ap 0 0 -beta];

B1 = [0;0;0;0;mi*gamma];

B2 = [0 0 0 0 0 0
    0 1/ms 0 0 0 0
    0 0 0 0 0 0
    -1 0 0 0 0 0
    0 0 0 0 0 0];

C = [1 0 0 0 0
    1 0 1 0 0
    -ks/ms -betas/ms 0 0 Ap/(mi*ms)];

D1 = [0;0;0];

D2 = [0 0 1 0 0 0
    0 0 0 1 0 0
    0 1/ms 0 0 1 0];

%% Initial Conditions: different from the equilibrium but not too much

Diseq = (0.1*(2*rand(1,1)-1));

x_init = [(- g*ms/ks + l0s + Lift/ks) - Diseq
    Diseq
    (Lift/kt + l0t - g*(mu+ms)/kt) - Diseq
    Diseq
    Diseq];

x_tilde_init = x_init - x0;

%% MATRICES OF e - e(1) = x(1) + vi(1) - r_l
Ce = [1 0 0 0 0];  

De1 = [0];

De2 = [0 0 1 0 0 -1];

%[V,Vn,J] = JCF(A);

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
    %Stabilizabilty Check
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
        eps5max = 1e6; % actuator force (do not pen)
        eps6max = 0.001; % integral of the position error PEN
        eps7max = 1e4; % sprung mass acceleration
        eps8max = 1e4; % sprung mass height 

     case 2
        eps1max = 0.001; % susp. deflection PEN
        eps2max = 0.001; % susp. speed PEN
        eps3max = 1e5; % tire deflection 
        eps4max = 1e5; % tire speed
        eps5max = 1e6; % actuator force (do not pen)
        eps6max = 0.001; % integral of the position error PEN
        eps7max = 1e4; % sprung mass acceleration
        eps8max = 0.001; % sprung mass height PEN
end

Q = inv(8*diag([eps1max^2,eps2max^2,eps3max^2,eps4max^2,eps5max^2, ...
    eps6max^2,eps7max^2,eps8max^2]));

umax = 0.01; % source: Road Adaptive....
R = inv(umax^2);
barR = R+D1eps.'*Q*D1eps;
alpha = 0; %we already have RE(lambda)>0. When alpha>6.6 Am,Ceps fully observable, otherwise is stabilisable
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
    imR = orth(R);              %pg.29 PDF toni
    imRorth = null(R.');         %ker(R)
    invTR = [imR imRorth];      %T_R^-1
    barA = inv(invTR)*Am*invTR;  
    A22 = barA(rankR+1:end, rankR+1:end);  %non-reachable part of barA
    eA22 = eig(A22);
    %Stabilizabilty Check
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
    %DECTABILITY CHECK
    A11 = barA(1:(length(Ceps)-rankO), 1:(length(Ceps)-rankO));  %non-observable part of barA
    eA11= eig(A11);
    %Stabilizabilty Check
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
XOinit = x_tilde_init;