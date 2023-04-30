function [d,zr] = disturbance(t,scenario)

%% BUMP MODEL
% scenario; % 0 = city; 1 = race; 2 = off roads

switch scenario
    case 0
        ds_dt = 30/3.6; % [m/s] vechile speed
        H = 0.10; % [m] bump height
        B = 1; % [m] bump base
    case 1
        ds_dt = 120/3.6; % [m/s] vechile speed
        H = 0.05; % [m] bump height
        B = 0.5; % [m] bump base
    case 2
        ds_dt = 3/3.6; % [m/s] vechile speed
        H = 0.5; % [m] bump height
        B = 2; % [m] bump base
    otherwise
        ds_dt = 0; % [m/s] vehicle speed
        H = 0.05; % [m] bump height
        B = 1; % [m] bump base
end

d2s_dt2 = 0; % [m/s^2] vechile acceleration
T = 5; % [s] bump time
L = ds_dt*T; % [m] bump location

c0 = 0;
c1 = 0;
c2 = 0;
c3 = 0;

x = [B^4 B^5 B^6 B^7 B^8
    (B/2)^4  (B/2)^5 (B/2)^6 (B/2)^7 (B/2)^8
    4*B^3 5*B^4 6*B^5 7*B^6 8*B^7
    12*B^2 20*B^3 30*B^4 42*B^5 56*B^6
    24*B 60*B^2 120*B^3 210*B^4 336*B^5]\[0; H; 0; 0; 0];
%c1 = x(1);
c4 = x(1);
c5 = x(2);
c6 = x(3);
c7 = x(4);
c8 = x(5);

s = t*ds_dt-L; % [m] travelled space
if s >=0 && s <= B
    zr = c0+c1*s+c2*s.^2+c3*s.^3+c4*s.^4+c5*s.^5+c6*s.^6+c7*s.^7+c8*s.^8; % [m] bump profile
    dzr_ds =c1+2*c2*s+3*c3*s.^2+4*c4*s.^3+5*c5*s.^4+6*c6*s.^5+7*c7*s.^6+8*c8*s.^7; % [-] bumb slope
    d2zr_ds2 = 2*c2+6*c3*s+12*c4*s.^2+20*c5*s.^3+30*c6*s.^4+42*c7*s.^5+56*c8*s.^6; % [1/m] bump curvature
    %d3zr_ds3 = 6*c3+24*c4*s+60*c5*s^2+120*c6*s^3+210*c7*s^4+336*c8*s^5; % [1/m] bump curvature
    %dot_zr = dzr_ds*ds_dt; % [m/s] bump speed
    ddot_zr = d2zr_ds2*(ds_dt).^2+dzr_ds.*d2s_dt2; % [m/s^2] bump acceleration
    d = ddot_zr;
else
    d = 0;
    zr = 0;
end





