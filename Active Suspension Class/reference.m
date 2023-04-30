function r = reference(t,scenario)

ks = 16000;
ms = 250;
ells = 0.5;
g = 9.81;

r0 = ells-g*ms/ks;
%scenario = 0; % 0 = city; 1 = race; 2 = off roads

switch scenario
    case 0
        Delta_r = 0;
    case 1
        Delta_r = -0.05;
    case 2
        Delta_r = 0.05;
    otherwise
        Delta_r = 0;
end


if t < 1
    r = r0;
else
    r = r0+Delta_r;
end

end