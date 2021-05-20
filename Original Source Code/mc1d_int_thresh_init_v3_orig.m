%Motor-clutch initialization

close all
clear all
clc

set(0,'DefaultFigureWindowStyle','docked')

nm = 150; %motors
Fm = 2; %stall force, pN
vu = 120; %unloaded velocity, nm/s
kaddb = 1; %basal cluch add-rate, s^-1
ncmax = 750; %max clutch number
Ft = 16; %reinforcement threshold, pN
nc = 150; %clutches
kon = 0.3; %on-rate for clutches, s^-1
koff = 0.1; %basal off-rate for clutches, s^-1
Fb = 2; %bond break force, pN
Fact = 2000; %initial F-actin, subunits
Gact = 2000; %initial G-actin, subunits
vp = 200; %leading edge polymerization velocity, nm/s

Kc = 0.8; %clutch stiffness, pN/nm
tmax = 10000; %max simulation time, s

Ksub = logspace(-1,2,20); %substrate stiffness, pN/nm

simdate = date;
mc_params_info = ['MC1D_retest_nm_' num2str(nm) '_nc_' num2str(nc)];
thresh_info = ['_Fthresh_' num2str(Ft) '_kaddb_' num2str(kaddb) ...
    '_ncmax_' num2str(ncmax) '_'];
time_info = ['simtime_' num2str(tmax)];

filename = [simdate mc_params_info thresh_info time_info '.mat'];

%Initialize
traction_force = [];
retrograde_flow = [];
clutch_num = [];
eng_clutch = [];
spread_ar = [];
Fperclutch = [];
Fperengclutch = [];

for ii = 1:size(Ksub,2)
    
    [Ftrac,vflow,nclutch,neng,cellar] = mc1d_int_thresh_func_v3_orig(tmax,...
        nm,Fm,vu,kaddb,ncmax,nc,kon,koff,Fb,Ft,Kc,Fact,Gact,vp,Ksub(ii));
    
    traction_force(ii) = Ftrac;
    retrograde_flow(ii) = vflow;
    clutch_num(ii) = nclutch;
    eng_clutch(ii) = neng;
    spread_ar(ii) = cellar;
    
    %Average forces
    Fperclutch(ii) = Ftrac/nclutch;
    Fperengclutch(ii) = Ftrac/neng;
    
    Ksub(ii)
    
end

Ksub = Ksub';
traction_force = traction_force';
retrograde_flow = retrograde_flow';
Fperclutch = Fperclutch';
Fperengclutch = Fperengclutch';
clutch_num = clutch_num';
eng_clutch = eng_clutch';
spread_ar = spread_ar';

figure()
subplot(3,3,1)
semilogx(Ksub,traction_force)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Traction Force (pN)')
axis([min(Ksub) max(Ksub) 0 (nm*Fm+5)])
subplot(3,3,2)
semilogx(Ksub,retrograde_flow)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Retrograde Flow (nm/s)')
axis([min(Ksub) max(Ksub) 0 vu])
subplot(3,3,3)
semilogx(Ksub,clutch_num)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Clutch Ensemble Size')
axis([min(Ksub) max(Ksub) 0 (max(clutch_num)+5)])
subplot(3,3,4)
semilogx(Ksub,eng_clutch)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Engaged Clutches')
axis([min(Ksub) max(Ksub) 0 (max(eng_clutch)+5)])
subplot(3,3,5)
semilogx(Ksub,spread_ar)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Cell Spread Area (\mum^2)')
axis([min(Ksub) max(Ksub) 0 (max(spread_ar)+500)])
subplot(3,3,6)
semilogx(Ksub,Fperclutch)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Force Per Clutch (pN)')
axis([min(Ksub) max(Ksub) 0 (max(Fperclutch+5))])
subplot(3,3,7)
semilogx(Ksub,Fperengclutch)
xlabel('Substrate Spring Constant (pN/nm)')
ylabel('Force Per Engaged Clutch (pN)')
axis([min(Ksub) max(Ksub) 0 (max(Fperengclutch)+5)])

save(filename)