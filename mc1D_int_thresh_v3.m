%1D Motor-clutch model v3
%August 25, 2015

%Uses direct Gillespie SSA, instead of indirect as in:
%Bangaser, Rosenfeld, & Odde, 2013

set(0,'DefaultFigureWindowStyle','docked')

%Integrin reinforcement - collaboration with UCSF
%1) Number of clutches influenced by kinetic equation that is scaled by the
%total available clutches in the cell
%   -kadd = kadd0*nc(F>Fthresh)*(ncmax-nc/ncmax)
%2) Reinforcement occurs when force is greater that >10pN on individual
%clutches (i.e. I is the number of clutches above this force threshhold and
%impacts the number of clutches added).

clear
clc
close all

%%
%Parameters

%Motors
n_motors = 150;
Force_motor = 2; %pN, stall force
velocity_unloaded = 120; %nm/s, unloaded velocity

%Clutches
k_add_base = 1; %Basal clutch addition rate (s^-1)
Force_thresh = 10; %pN, threshold force for integrin reinforcement
n_clutch_max = 450; %
n_clutch = 75;
k_on = 0.3; %s-1, pseudo-first order clutch on-rate
k_off = 0.1; %s-1, basal first-order off-rate
Force_bond = 2; %pN, characteristic break force
stiffness_clutch = 0.8; %pN/nm, clutch stiffness

%%
%Length Parameters
nucleus = 5000; %nm, nuclear radius
F_actin(1) = 2000; %Length in subunits
Length_actin = 4; %nm, dimer length
Length_cell(1) = nucleus + (F_actin(1)*Length_actin);
velocity_polymer = 200; %nm/s, leading edge polymerization velocity
G_actin(1) = 2000; %subunits, G-actin pool
total_actin(1) = F_actin(1)+G_actin(1);
% clutchpos(1) = Length_cell; %Position of clutches relative to cell edge

%Substrate
stiffness_substrate = 1; %pN/nm, substrate stiffness

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
max_time = 5000;

% clutch_rates = [];
traction_force = [];
retrograde_flow = [];
engaged_clutches = [];
n_thresh = [];
clutch_add_rate = [];
clutch_forces_nz_max = [];
cell_area = [];
clutchsize = n_clutch;
cycles = 0;

%%
while simtime < max_time
    
    %Calculate off-rate for engaged clutches based on clutch deformations
    clutch_Forces = stiffness_clutch.*(x_clutch-x_substrate);
    norm_clutch_Forces = clutch_Forces./Force_bond;
    clutch_unbinding = k_off.*exp(norm_clutch_Forces);
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
        clutch_rates(n_clutch) = k_off;
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
    clutchsize(event) = n_clutch;
%     cell_area(event) = pi*(Length_cell(event)^2);
%     clutchpos(event) = Length_cell(event)-x_substrate;
    n_thresh(event) = size(clutch_Forces_threshold,1);
    t_step(event) = simtime(event)-simtime(event-1);
    
    %Count cycles
    if x_substrate == 0
        cycles = cycles+1;
    end
    
    %Nonzero mean of clutch forces
    clutch_Forces_nz_max(event) = max(clutch_Forces);
    
end

%Numerical outputs
TF = sum((traction_force.*t_step)/simtime(end))
RF = sum((retrograde_flow.*t_step)/simtime(end))
cycle_time = simtime(end)/cycles
Force_bins = 0:20;
max_Force_engaged_clutches = hist(clutch_Forces_nz_max,Force_bins);
norm_max_Force_engaged_clutches = max_Force_engaged_clutches./event;
cell_area = (pi*Length_cell.^2).*1e-6;

figure()
%Substrate position
subplot(4,2,1)
plot(simtime,(-substrate_position))
xlabel('Time (s)')
ylabel('X-position (nm)')
axis([0 simtime(end) (-max(substrate_position)-100) 0])
%Number of clutches in ensemble
subplot(4,2,2)
plot(simtime,clutchsize)
xlabel('Time (s)')
ylabel('Total Clutches')
axis([0 simtime(end) 0 max(clutchsize)+5])
%Number of engaged clutches
subplot(4,2,3)
plot(simtime,engaged_clutches)
xlabel('Time (s)')
ylabel('Engaged Clutches')
axis([0 simtime(end) 0 max(clutchsize)+5])
%Thresholded clutches
subplot(4,2,4)
plot(simtime,n_thresh)
xlabel('Time (s)')
ylabel('Thresholded Clutches')
axis([0 simtime(end) 0 max(clutchsize)+5])
%Cell Area
subplot(4,2,5)
plot(simtime,cell_area)
xlabel('Time (s)')
ylabel('Cell Spread Area (\mum^2)')
axis([0 simtime(end) 0 (max(cell_area)+500)])
%F-actin and G-actin
subplot(4,2,6)
plot(simtime,F_actin,'b',simtime,G_actin,'r',simtime,total_actin,'k')
xlabel('Time (s)')
ylabel('Actin (subunits)')
axis([0 simtime(end) 0 (max(total_actin)+100)])
%Integrin max force sensor data
subplot(4,2,7)
bar(Force_bins,norm_max_Force_engaged_clutches)
xlabel('Force (pN)')
ylabel('Frequency')
axis([0 max(Force_bins) 0 1])