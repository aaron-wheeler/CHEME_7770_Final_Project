%% Cell Spread Area Plot - Final Project - CHEME 7770 - UPDATED: 20210522

% Note: All variables were imported from corresponding environments
% manually and script was adapted from Mekhdjian et al. 2017

figure(1)
plot(simtime_k1rup1,cell_area_k1rup1)
hold on
plot(simtime_k1rup10,cell_area_k1rup10)
plot(simtime_k1rup150,cell_area_k1rup150)
hold off
title('Cell Spread Area for a Range of k_1_,_r_u_p Values')
xlabel('Time (s)')
ylabel('Cell Spread Area (\mum^2)')
legend('k_1_,_r_u_p = 1','k_1_,_r_u_p = 10','k_1_,_r_u_p = 150','Location','southeast')
axis([0 5000 0 1000])