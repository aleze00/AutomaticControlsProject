Res = load("OffRoad_Obs_bump.mat");   %salvato con save("ciccio","out")

%% Plot EPS_3
figure();
plot(Res.out.y_CL_NL.Time, Res.out.y_CL_NL.Data(3,:) , 'g', 'DisplayName','Acceleration [m/s^2]');
hold on
plot(Res.out.y_CL_NL.Time, Res.out.y_CL_NL.Data(2,:) , 'b', 'DisplayName','Suspension Length [m]');
plot(Res.out.y_CL_NL.Time, Res.out.y_CL_NL.Data(1,:), 'r', 'DisplayName','Suspension Deflection [m]');
legend();
title(strcat("Outputs using Bump without Observer"));
xlabel("simulation time [s]")
ylabel("")
hold off;
