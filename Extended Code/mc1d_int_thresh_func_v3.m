function [TF,RF,nc,neng,cellar] = mc1d_int_thresh_func_v3(max_time,...
    n_motors,Force_motor,velocity_unloaded,k_add_base,n_clutch_max,...
    n_clutch,k_on,k_off,Force_bond,Force_thresh,stiffness_clutch,...
    F_actin_init,G_actin_init,velocity_polymer,stiffness_substrate)

%1D Motor-clutch model v3
%Aug 13, 2015

%Uses direct Gillespie SSA, instead of indirect as in:
%Bangaser, Rosenfeld, & Odde, 2013

%Integrin reinforcement - collaboration with UCSF
%1) Number of clutches determined by kinetic equation k = k0(Ncmax-Nc),
%where k0 is determined by the force per bond
%2) Reinforcement occurs when force is greater that >2pN on individual
%clutches (i.e. I is the number of clutches above this force threshhold and
%impacts the number of clutches added).

%Fixed 9/9/15 - error where spread area had a minimum at ~1pN

% clear
% clc
% close all

%%
% %Parameters
% 
% %Motors
% n_motors = 50;
% Force_motor = 2; %pN, stall force
% velocity_unloaded = 120; %nm/s, unloaded velocity
% 
% %Clutches
% k_in_base = 1; %Basal clutch addition rate (s^-1)
% k_out = 0.02; %Basal clutch removal rate (s^-1)
% thresh_gain = 0.003; %Gain factor for integrin addition
% n_clutch = round(k_in_base/k_out);
% k_on = 0.3; %s-1, pseudo-first order clutch on-rate
% k_off = 0.1; %s-1, basal first-order off-rate
% Force_bond = 2; %pN, characteristic break force
% Force_thresh = 1.5; %pN, threshold force for integrin reinforcement
% stiffness_clutch = 0.8; %pN/nm, clutch stiffness
% 
% %Substrate
% stiffness_substrate = 0.1; %pN/nm, substrate stiffness

%%
%Initialization
x_substrate = 0;
x_clutch = zeros(n_clutch,1);
Force_substrate = stiffness_substrate*x_substrate;
clutch_state = zeros(n_clutch,1); %0 unbound, 1 bound, start all unbound
clutch_Forces = stiffness_clutch.*(x_clutch-x_substrate);
n_engaged = sum(clutch_state);
event = 1;
simtime(1) = 0;
% max_time = 5000;

%%
%Length Parameters
nucleus = 5000; %Nuclear radius
F_actin(1) = F_actin_init; %Length in subunits
Length_actin = 4; %nm, dimer length
Length_cell(1) = nucleus+F_actin*Length_actin;
G_actin(1) = G_actin_init; %subunits, G-actin pool
total_actin(1) = F_actin(1)+G_actin(1);
% cell_area(1) = pi*Length_cell(1)^2*1e-6; %um^2, "area" of a circular cell
% clutchpos(1) = Length_cell; %Position of clutches relative to cell edge

clutch_rates = [];
traction_force = [];
retrograde_flow = [];
engaged_clutches = [];
n_thresh = [];
clutch_add_rate = [];
clutchsize = n_clutch;
cycles = 0;
delE21 = (4.8*1000*(4.114/9.83E-22)) / 6.022E23; %kcal/mol * cal/kcal * pN*nm/cal / mol^-1 = pN*nm Difference in energy level between folded and unfolded state in solution under constant T,P
kbT = 4.114*(310/298); %pN*nm    boltzman constant * 310K
phi = exp(delE21/kbT); % Equilibrium Occupancy Parameter ~ Occ @ zero force
k1rup = 10; % dissociation rate of off rate at zero force  [s^-1]
Fu = 30; % Force of 1/2 occupancy of the unfolded and folded states under AFM [pN]
f12 = Fu/(delE21/kbT); % scaling force for low force pathway [pN]

%%
while simtime < max_time
   
    %Calculate off-rate for engaged clutches based on clutch deformations
    clutch_Forces = stiffness_clutch.*(x_clutch-x_substrate);
    norm_clutch_Forces = clutch_Forces./Force_bond;
    clutch_unbinding = (phi*k1rup + exp(clutch_Forces./f12) .* k_off.*exp(norm_clutch_Forces)) ./ (phi + exp(clutch_Forces./f12));
    clutch_unbinding = clutch_unbinding.*(clutch_state == 1);
    clutch_binding = k_on.*(clutch_state == 0);
    clutch_rates = clutch_unbinding+clutch_binding;
    
    %Calculate add/loss rates for individual clutches
    clutch_Forces_threshold = clutch_Forces(clutch_Forces>Force_thresh);
    k_add = k_add_base*size(clutch_Forces_threshold,1)*...
        ((n_clutch_max-n_clutch)/n_clutch_max);
    clutch_rates(n_clutch+1) = k_add; %adding clutch rate
    
    %Determine event time based on clutch on- and off-rates
    URN1 = rand;
    event = event+1;
    event_time = -log(URN1)/sum(clutch_rates);
    simtime(event) = simtime(event-1)+event_time;
    
    %Normalized rates and event selection
    normalized_rate = cumsum(clutch_rates)./sum(clutch_rates);
    URN2 = rand;
    selected_event = find((normalized_rate>URN2),1,'first');
     
    %Add/loss events
    if selected_event == n_clutch+1; %Addition
        n_clutch = n_clutch+1;
        x_clutch(n_clutch) = 0;
        clutch_state(n_clutch) = 0;
        clutch_Forces(n_clutch) = 0;
        clutch_rates(n_clutch) = (phi*k1rup + k_off) ./ (phi + 1);
    %Execute the event (binding or unbinding of selected clutch)
    elseif clutch_state(selected_event) == 0 %unbound
        clutch_state(selected_event) = 1;
        x_clutch(selected_event) = x_substrate;
    elseif clutch_state(selected_event) == 1 %bound
        clutch_state(selected_event) = 0;
    end
    
    %Calculate F-actin position change by retrograde flow
    Force_filament = Force_substrate/(Force_motor*n_motors);
    filament_velocity = velocity_unloaded*(1-Force_filament);
    
    %Advance engaged clutch positions
    clutch_advance = filament_velocity*event_time;
    x_clutch = x_clutch+(clutch_advance.*(clutch_state == 1));
    
    %Advance cell edge position and cell length
    velocity_edge = velocity_polymer*((G_actin(event-1)/(total_actin(event-1))));
    edge_advance = velocity_edge*event_time;
%     Length_cell(event) = Length_cell(event-1)+Length_actin*round((edge_advance-clutch_advance)/Length_actin);
    
    %Calculate change in F-actin and G-actin
    G_actin_gain = round(clutch_advance/Length_actin);
    G_actin_loss = round(edge_advance/Length_actin);
    F_actin(event) = F_actin(event-1)-G_actin_gain+G_actin_loss;
    G_actin(event) = G_actin(event-1)+G_actin_gain-G_actin_loss;
    total_actin(event) = F_actin(event)+G_actin(event);
    Length_cell(event) = nucleus + (F_actin(event)*Length_actin);
    
    %Calculate substrate position by force balance
    traction_force(event) = Force_substrate;
    n_engaged = sum(clutch_state);
    sum_clutch_Forces = stiffness_clutch.*sum(x_clutch(clutch_state==1));
    x_substrate = sum_clutch_Forces./(stiffness_substrate+(n_engaged*stiffness_clutch));
    Force_substrate = stiffness_substrate*x_substrate;
    x_clutch(clutch_state==0) = x_substrate;
    
    %Calculate/tabulate outputs
    substrate_position(event) = x_substrate;
    engaged_clutches(event) = n_engaged;
    retrograde_flow(event) = filament_velocity;
%     cell_area(event) = pi*(Length_cell(event)^2)*1e-6;
    clutchsize(event) = n_clutch;
    n_thresh(event) = size(clutch_Forces_threshold,1);
    t_step(event) = simtime(event)-simtime(event-1);
    
    %Count cycles
    if x_substrate == 0
        cycles = cycles+1;
    end
    
end

cell_area = (pi*Length_cell.^2).*1e-6;

% figure()
% %Substrate position
% subplot(3,2,1)
% plot(simtime,(-substrate_position))
% title('Substrate Position')
% xlabel('Time (s)')
% ylabel('X-position (nm)')
% axis([0 simtime(end) (-max(substrate_position)-100) 0])
% %Number of clutches in ensemble
% subplot(3,2,2)
% plot(simtime,clutchsize)
% title('Clutch Number')
% xlabel('Time (s)')
% ylabel('Total Clutches')
% axis([0 simtime(end) 0 max(clutchsize)+5])
% %Number of engaged clutches
% subplot(3,2,3)
% plot(simtime,engaged_clutches)
% title('Number of Engaged Clutches')
% xlabel('Time (s)')
% ylabel('Engaged Clutches')
% axis([0 simtime(end) 0 max(clutchsize)+5])
% %Thresholded clutches
% subplot(3,2,4)
% plot(simtime,n_thresh)
% title('Clutches Above Threshold Force')
% xlabel('Time (s)')
% ylabel('Thresholded Clutches')
% axis([0 simtime(end) 0 max(clutchsize)])
% %Clutch addition rate
% subplot(3,2,5)
% plot(simtime,clutch_add_rate)
% title('Clutch Ensemble Addition Rate')
% xlabel('Time (s)')
% ylabel('Clutch Addition Rate (1/s)')
% axis([0 simtime(end) (1.1*min(clutch_add_rate)) (1.1*max(clutch_add_rate))])

%Numerical outputs
TF = sum((traction_force.*t_step)/simtime(end));
RF = sum((retrograde_flow.*t_step)/simtime(end));
nc = sum((clutchsize.*t_step)/simtime(end));
neng = sum((engaged_clutches.*t_step)/simtime(end));
cellar = sum((cell_area.*t_step)/simtime(end));
% cycle_time = simtime(end)/cycles
end