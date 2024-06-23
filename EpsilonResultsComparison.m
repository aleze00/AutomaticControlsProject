%% Load Results of Simulation
res0_NoObs = load("res0race.mat");
res1_NoObs = load("res8race1.mat");
res2_NoObs = load("res8race2.mat");
res3_NoObs = load("Final3_Comfort.mat");
res4_NoObs = load("Final4_Comfort.mat");
res5_NoObs = load("Final5_Comfort.mat");
res6_NoObs = load("Final6_Comfort.mat");
res7_NoObs = load("Final7_Comfort.mat");
res8_NoObs = load("Final8_Comfort.mat");

%% Plot Comparison
for i = 1 : 8 
    plotComparison1(res0_NoObs,res1_NoObs,res2_NoObs, i);
    % plotComparison(res0_NoObs,res1_NoObs,res2_NoObs,res3_NoObs,res4_NoObs,res5_NoObs,res6_NoObs,res7_NoObs,res8_NoObs, i);
end

function plotComparison1(res0,res1,res2,index_eps)

figure();
plot(res0.out.eps.Time, res0.out.eps.Data(:,index_eps), 'r', 'DisplayName','Res0');
hold on
plot(res1.out.eps.Time, res1.out.eps.Data(:,index_eps) , 'g', 'DisplayName','Res1');
plot(res2.out.eps.Time, res2.out.eps.Data(:,index_eps), 'b', 'DisplayName','Res2');
legend();
title(strcat("Eps_", string(index_eps)));
hold off;

end

function plotComparison(res0,res1,res2,res3,res4,res5,res6,res7,res8,index_eps)

figure();
plot(res0.out.eps.Time, res0.out.eps.Data(:,index_eps), 'r', 'DisplayName','Res0');
hold on
plot(res1.out.eps.Time, res1.out.eps.Data(:,index_eps) , 'g', 'DisplayName','Res1');
plot(res2.out.eps.Time, res2.out.eps.Data(:,index_eps), 'b', 'DisplayName','Res2');
plot(res3.out.eps.Time, res3.out.eps.Data(:,index_eps) , 'c', 'DisplayName','Res3');
plot(res4.out.eps.Time, res4.out.eps.Data(:,index_eps), 'm', 'DisplayName','Res4');
plot(res5.out.eps.Time, res5.out.eps.Data(:,index_eps) , 'Color', '#4DBEEE', 'DisplayName','Res5');
plot(res6.out.eps.Time, res6.out.eps.Data(:,index_eps), 'k', 'DisplayName','Res6');
plot(res7.out.eps.Time, res7.out.eps.Data(:,index_eps) , 'Color', '#EDB120', 'DisplayName','Res7');
plot(res8.out.eps.Time, res8.out.eps.Data(:,index_eps), 'Color', '#7E2F8E', 'DisplayName','Res8');
legend();
title(strcat("Eps_", string(index_eps)));
hold off;

end