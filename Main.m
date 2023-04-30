clc
close all
clear all

%% INITIAL CONDITIONS
x0 = [0.01
    0
    0.01
    0];

%% MODEL INVESTIGATION

kt = 160000; % [N/m] tire stiffness

mu = 25; % [kg] unsprung mass
ms = 300; % [kg] sprung mass
barm = 1/mu + 1/ms; % [1/kg]

betas = 1000; % [N/(m/s)]  damping coeff.
ks = 16000; % [N/m] suspension stiffness

A = [ 0 1  0 0
     -ks*barm-kt/mu -betas*barm kt/mu 0
     0 0 0 1
     -ks/ms -betas/ms 0 0];


%% EQUILIBRIUM POINT
ls = 0.5; % [m] suspension unloaded length
lt = 0.05; % [m] tire unloaded length
g = 9.83; % [m/s^2] gravity acceleration

B3 = [0 0 0
      ks*barm -kt/mu 0
      0 0 0
      ks/ms 0 -1];

xeq = -inv(A)*B3*[ls; lt; g];

%% REACHABILITY CHECK

B1 = [0; barm; 0; 1/ms];
R = ctrb(A,B1);
rankR = rank(R);
if rankR == length(A)
    disp('FULLY REACHABLE')
else
    disp('NOT FULLY REACHABLE')
    imR = orth(R);
    imRorth = ker(R.');
    invTR = [imR imRorth];
    barA = inv(invTR)*A*invTR;
    A22 = barA(rankR+1:end, rankR+1:end);
    eA22 = eig(A22);
    for i = 1:length(eA22)
        if real(eA22(i)) >= 0
            disp('NOT STABILISABLE')
        end
    end
end

%% OBSERVABILITY CHECK

C1 = [1 0 0 0
      -ks/ms -betas/ms 0 0];
  
O = obsv(A,C1);
rankO = rank(O);

if rankO == length(A)
    disp('FULLY OBSERVABLE')
else
    disp('NOT FULLY OBSERVABLE')
end