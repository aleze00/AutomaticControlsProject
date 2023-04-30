function [dot_x,y,e] = F(x,u,w)

%% VARAIBLES

% STATES
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

% DISTURBANCE
d = w(1);

% NOISE
nu = w(2:3);

% REFERENCE
r = w(4);

%% ACTUAL PLANT PARAMETERS

ks = 16000;
bs = 1000;
mu = 45;
ms = 250;
ells = 0.5;
barm = (mu+ms)/ms;
g = 9.81;

%% STATE DYNAMICS

dot_x = [x2
    barm/mu*(-ks*(x1-ells)-bs*x2+u)-1/mu*ft(x3)
    x4
    -g+ks/mu*(x1-ells)+bs/mu*x2-u/mu+1/mu*ft(x3)-d];

%% MEASUREMENTS

y = [x1;
    -g-ks/ms*(x1-ells)-bs/ms*x2+u/ms+g] + nu;

%% CONTROLLED VARIABLES

e = x1-r;

end

%% SUPPORT FUNCTIONS

function f = ft(s)
ell = 0.10; 
kt = 10*16000;
if s > ell
    f = 0;
else
    f = -kt*(s-ell);
end
end