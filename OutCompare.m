Res = load("Race_Obs_Road.mat");   %salvato con save("ciccio","out")

%% Plot EPS_3
figure();
plot(Res.out.y_CL_NL.Time, Res.out.y_CL_NL.Data(1,:), 'r', 'DisplayName','Suspension Deflection [m]');
hold on
plot(Res.out.y_CL_NL.Time, Res.out.y_CL_NL.Data(2,:) , 'b', 'DisplayName','Suspension Length [m]');
plot(Res.out.y_CL_NL.Time, Res.out.y_CL_NL.Data(3,:) , 'g', 'DisplayName','Acceleration [m/s^2]');
legend();
title(strcat("Outputs using Road Disturbance with Observer"));
xlabel("simulation time [s]")
ylabel("")
hold off;
