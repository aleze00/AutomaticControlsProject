%% Load Results of Simulation
res0_NoObs = load("Res_NO_OBS_DM_0_Comfort.mat");
res1_NoObs = load("Res_NO_OBS_DM_1_OffRoad.mat");
res2_NoObs = load("Res_NO_OBS_DM_2_Race.mat");
%% Plot Comparison
for i = 1 : 8 
    plotComparison(res0_NoObs,res1_NoObs,res2_NoObs, i);
end


function plotComparison(res0,res1,res2,index_eps)

figure();
plot(res1.out.eps.Time, res1.out.eps.Data(:,index_eps), 'g', 'DisplayName','OffRoad');
hold on
plot(res0.out.eps.Time, res0.out.eps.Data(:,index_eps) , 'r', 'DisplayName','Comfort');
plot(res2.out.eps.Time, res2.out.eps.Data(:,index_eps), 'b', 'DisplayName','Race');
legend();
title(strcat("Eps_", string(index_eps)));
hold off;

end