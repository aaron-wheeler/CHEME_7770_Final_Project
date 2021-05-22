%% Clutch Slip across Ksub Comparison Plots - Final Project - CHEME 7770 - UPDATED: 20210522

% Note: All variables were imported from corresponding environments
% manually and script was adapted from Mekhdjian et al. 2017
% Variable Renaming Commands used between environment imports:
% Ksub_cs = Ksub;
% traction_force_cs = traction_force;
% retrograde_flow_cs = retrograde_flow;
% clutch_num_cs = clutch_num;
% eng_clutch_cs = eng_clutch;
% spread_ar_cs = spread_ar;
% Fperclutch_cs = Fperclutch;
% Fperengclutch_cs = Fperengclutch;

figure()
subplot(3,3,1)
semilogx(Ksub,traction_force)
hold on
semilogx(Ksub_cs,traction_force_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Traction Force (pN)')
axis([min(Ksub) max(Ksub) 0 (nm*Fm+5)])
legend('Traditional Bond','Clutch-Slip Bond')

subplot(3,3,2)
semilogx(Ksub,retrograde_flow)
hold on
semilogx(Ksub_cs,retrograde_flow_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Retrograde Flow (nm/s)')
axis([min(Ksub) max(Ksub) 0 vu])
legend('Traditional Bond','Clutch-Slip Bond','Location','southeast')

subplot(3,3,3)
semilogx(Ksub,clutch_num)
hold on
semilogx(Ksub_cs,clutch_num_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Clutch Ensemble Size')
axis([min(Ksub) max(Ksub) 0 (max(clutch_num_cs)+10)])
legend('Traditional Bond','Clutch-Slip Bond','Location','southeast')

subplot(3,3,4)
semilogx(Ksub,eng_clutch)
hold on
semilogx(Ksub_cs,eng_clutch_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Engaged Clutches')
axis([min(Ksub) max(Ksub) 0 (max(eng_clutch)+5)])
legend('Traditional Bond','Clutch-Slip Bond')

subplot(3,3,5)
semilogx(Ksub,spread_ar)
hold on
semilogx(Ksub_cs,spread_ar_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Cell Spread Area (\mum^2)')
axis([min(Ksub) max(Ksub) 0 (max(spread_ar)+500)])
legend('Traditional Bond','Clutch-Slip Bond')

subplot(3,3,6)
semilogx(Ksub,Fperclutch)
hold on
semilogx(Ksub_cs,Fperclutch_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Force Per Clutch (pN)')
axis([min(Ksub) max(Ksub) 0 (max(Fperclutch+5))])
legend('Traditional Bond','Clutch-Slip Bond')

subplot(3,3,7)
semilogx(Ksub,Fperengclutch)
hold on
semilogx(Ksub_cs,Fperengclutch_cs)
hold off
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Force Per Engaged Clutch (pN)')
legend('Traditional Bond','Clutch-Slip Bond')
axis([min(Ksub) max(Ksub) 0 (max(Fperengclutch)+5)])