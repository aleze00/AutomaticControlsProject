clc 
clear 
close all

% (x1,x2,x3,x4):=(zs,˙vs,zu,˙vu)  

% Physical parameters
ms = 300;    % kg
mu = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

% State matrices
A = [ 0 1 0 0; 
    [-ks -bs ks bs]/ms ;
      0 0 0 1; 
      [ks bs -ks-kt -bs]/mu];

B = [ 0 0; 
    0 1e3/ms ; 
    0 0 ; 
    [kt -1e3]/mu];

C = [1 0 0 0; 
    1 0 -1 0; 
    A(2,:)];

D = [0 0; 
    0 0; 
    B(2,:)];

qcar = ss(A,B,C,D);

qcar.StateName = {'zs (m)';'vs (m/s)';...
          'zu (m)';'vu (m/s)'};
qcar.InputName = {'r';'fs'};
qcar.OutputName = {'zs';'SusLen';'AccS'};

bodemag(qcar({'zs','SusLen'},'r'),'b',qcar({'zs','SusLen'},'fs'),'r',{1 100});
legend('Road disturbance (r)','Actuator force (fs)','location','SouthWest')
title({'Gain from road dist (r) and actuator force (fs) ';
       'to body accel (AccS) and suspension travel (SusLen)'})

Wunc = makeweight(0.40,15,3);
unc = ultidyn('unc',[1 1],'SampleStateDim',5);
Act = ActNom*(1 + Wunc*unc);
Act.InputName = 'u';
Act.OutputName = 'fs';

