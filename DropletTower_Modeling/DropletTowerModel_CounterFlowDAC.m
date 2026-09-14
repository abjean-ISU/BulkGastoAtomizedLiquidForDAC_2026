clear all;close all;clc;

%DAC droplet tower model with water

%USER INPUTS---------------------------------------------------------------

%Import droplet size data OR Designate mean droplet diameter (mm) - 
% select relavent input
DropletSize = readmatrix(['C:\Users\anbarron\Box\Research\Projects\DAC\Droplet Tower Modeling\DropletDiameterData.csv']); %droplet diameter data (mm) 
%DropletSize = 1.6; %mean droplet diameter (mm)

%Input Tower Parameters 
towerHeight = 3; %m
towerDiameter = 12; %m 

%Input manifold parameters
d_orifice = 0.8/1000; %droplet manifold orifice diameter (m)
N_o = 1112203; %number of orifices in droplet manifold 

%Input operating conditions 
Q = 4.17E-01; %liquid flowrate (m3/s)
T = 20+273.15; %liquid temperature (K)
P = 1; %operating (headspace) pressure (atm) 

%Gas Composition 
DAB_CO2 = 1.67E-9; %diffusion coefficient of CO2 in water at correct Temp (m2/s)
DAB_O2 = 2.01E-9; %diffusion coefficient of O2
DAB_N2 = 2.00E-9; %diffusion coefficient of N2
MW_CO2 = 44.01; %CO2 gas molar mass (g/mol)
MW_O2 = 32; %O2 gas molar mass (g/mol)
MW_N2 = 28; %N2 gas molar mass (g/mol)
MW_Air = 28.97; %Air molar mass (g/mol)
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
ppm0_CO2 = 400; %initial CO2 conentration (ppm)
ppm0_O2 = 210000; %initial O2 conentration (ppm)
ppm0_N2 = 780000; %initial N2 conentration (ppm)
Ug_range = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0];  % m/s

c0_CO2 = 0; %initial CO2 dissolved gas concentration (g/L)
c0_O2 = 0; %initial O2 dissolved gas concentration (g/L)
c0_N2 = 0; %initial N2 dissolved gas concentration (g/L)

%Constant parameters 
g = 9.8; %acceleration due to gravity (m/s2)
rho_g = 1.184; %density of air (kg/m3)
rho_l = 997; %density of water (kg/m3) 
mu_g = 1.53E-5; %viscocity of air (m2/s)
s_tension = 72.8; %water surface tenstion (g/s2)

%END USER INPUTS ----------------------------------------------------------

Cg_t = P*MW_Air/(R*T)/1000; %total concentraiton of air (g/L)

H_CO2 = (1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T)))/55342*1000/MW_CO2; %Henry's Constant for CO2 (atm-L/g)
H_O2 = (1/(exp(-171.2542+8391.24/T+23.24323*log(T)+0*T)))/55342*1000/MW_O2; %Henry's Constant for O2 (atm-L/g)
H_N2 = (1/(exp(-181.587+8632.13/T+24.7981*log(T)+0*T)))/55342*1000/MW_N2; %Henry's Constant for N2 (atm-L/g)

%Determining weights and abscissae for droplet size distribution
if size(DropletSize)  == 1
    GQMOM_model_x = DropletSize/1000/2; %convert mean drop diameter to radius (m)
    GQMOM_model_w = 1; %using only 1 abscissae
else 
    %normalizing diameters against maximum diameter 
    diameter_max = max(DropletSize);
    diameter_norm = DropletSize./diameter_max;

    %computing the histogram and moments
    [hist_weights,edges] = histcounts(diameter_norm,20,'Normalization','probability');
    hist_abscissae = edges(1:end-1)+ diff(edges)/2;
    for i = 0:10
        calc = hist_weights.*(hist_abscissae.^i);
        hist_moments(i+1) = sum(calc);
    end

    %Using GQMOM with 4 abscissae to determine ODEs 
    [GQMOM_model_w,GQMOM_model_x,GQMOM_model_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),4);

    %converting x to radii in m
    GQMOM_model_x = GQMOM_model_x*diameter_max/2/1000; %m
end

M = length(GQMOM_model_x); %number of quadrature moments

%Computing initial system parameters
%height
xspan = 0:0.01:towerHeight; %location (m)
%velocity 
A_orifice = pi*(d_orifice/2)^2*N_o; %cross sectional orifice area (m2)
A_c = pi*(towerDiameter/2)^2; %column cross sectional area (m2)
v0 = Q/A_orifice; %initial droplet velocity (m/s)
%time
t0 = 0; %inital time (s)

%intial conditions vector
y0 = [repmat([v0,c0_CO2,c0_O2,c0_N2,t0],1,M),ppm0_CO2,ppm0_O2,ppm0_N2];

%BVP solver options
options = bvpset('RelTol', 1e-4, 'AbsTol', 1e-6, 'NMax', 2000);

flow_configs = {'coflow','counterflow'};

%Preallocate results - 3D: [flow config, Ug, x]
ppm_CO2_out_all  = zeros(2, length(Ug_range), length(xspan));
ppm_O2_out_all   = zeros(2, length(Ug_range), length(xspan));
ppm_N2_out_all   = zeros(2, length(Ug_range), length(xspan));
c_CO2_out_all    = zeros(2, length(Ug_range), length(xspan));
c_O2_out_all     = zeros(2, length(Ug_range), length(xspan));
c_N2_out_all     = zeros(2, length(Ug_range), length(xspan));
v_out_all        = zeros(2, length(Ug_range), length(xspan));
t_out_all        = zeros(2, length(Ug_range), length(xspan));
capture_eff_all  = zeros(2, length(Ug_range));
CO2_purity       = zeros(2, length(Ug_range));
kLa_L_CO2_all   = zeros(2, length(Ug_range));
kLa_T_CO2_all   = zeros(2, length(Ug_range));
kLa_L_O2_all    = zeros(2, length(Ug_range));
kLa_T_O2_all    = zeros(2, length(Ug_range));
kLa_L_N2_all    = zeros(2, length(Ug_range));
kLa_T_N2_all    = zeros(2, length(Ug_range));
holdup_all       = zeros(2, length(Ug_range));
CO2_captured_all = zeros(2, length(Ug_range));

for f = 1:2
    for j = 1:length(Ug_range)
        
        Ug = Ug_range(j); 
        lambda = 10^6/(Ug*Cg_t); %conversion from mdot to ppm

        constants = [g, rho_g, rho_l, mu_g, DAB_CO2, DAB_O2, DAB_N2, s_tension, P, H_CO2, H_O2, H_N2, A_c, Q, lambda]; 
        
        if strcmp(flow_configs{f}, 'coflow')
            %co-current - IVP solved with ode45
            [~, y_out] = ode45(@(x,y) ODEsys_Angelo_CoFlowcoupled(x, y,GQMOM_model_x, GQMOM_model_w, constants), xspan, y0);
        
        else 
        %counter-flow - BVP solved with bvp4c
            if j == 1
                solinit = bvpinit(xspan, y0);
            else
                solinit = bvpinit(xspan, @(x) interp1(xspan, y_out, x)');
            end
            sol   = bvp4c(@(x,y) ODEsys_Angelo_CounterFlowcoupled(x, y, GQMOM_model_x, GQMOM_model_w, constants), @(ya,yb) BC_coupled(ya, yb, v0, ppm0_CO2, ppm0_O2, ppm0_N2, M),solinit, options);
            y_out = deval(sol, xspan)';
        end    
        
        for i = 1:M
            idx      = (i-1)*5 + 1;
            v_ri(:,i) = y_out(:,idx);
            c_CO2_ri(:,i) = y_out(:,idx+1);
            c_O2_ri(:,i) = y_out(:,idx+2);
            c_N2_ri(:,i) = y_out(:,idx+3);
            t_ri(:,i) = y_out(:,idx+4);
        end

        %QMOM weighted averages
        v_x   = sum(v_ri .* GQMOM_model_w', 2);
        c_CO2_x   = sum(c_CO2_ri .* GQMOM_model_w', 2);
        c_O2_x   = sum(c_O2_ri .* GQMOM_model_w', 2);
        c_N2_x   = sum(c_N2_ri .* GQMOM_model_w', 2);
        t_x   = sum(t_ri .* GQMOM_model_w', 2);
        ppm_CO2_x = y_out(:,end-2);
        ppm_O2_x = y_out(:,end-1);
        ppm_N2_x = y_out(:,end);
    
         %store profiles
        ppm_CO2_out_all(f,j,:) = ppm_CO2_x;
        ppm_O2_out_all(f,j,:) = ppm_O2_x;
        ppm_N2_out_all(f,j,:) = ppm_N2_x;
        c_CO2_out_all(f,j,:)   = c_CO2_x;
        c_O2_out_all(f,j,:)   = c_O2_x;
        c_N2_out_all(f,j,:)   = c_N2_x;
        v_out_all(f,j,:)   = v_x;
        t_out_all(f,j,:)   = t_x;
        
        %system parameters 
        tower_v_tot = pi*(towerDiameter/2)^2*towerHeight; %tower volume (m3)
        
        %droplet parameters
        droplet_r = sum(GQMOM_model_x.*GQMOM_model_w); %average droplet radius (m)
        droplet_v = 4/3*pi*droplet_r^3; %average droplet volume (m3)
        
                %Liquid holdup parameters
        droplet_v_tot_m3 = Q*max(t_x); %total droplet (LH) volume (m3)
        droplet_v_tot_L = 1000*droplet_v_tot_m3; %total droplet (LH) volume (L)
        N_tot = droplet_v_tot_m3/droplet_v; %total number of droplets
        droplet_density_tot = N_tot/tower_v_tot; %total droplet number density (droplets/m3)
        liquid_holdup_tot = droplet_v_tot_m3/tower_v_tot; %total liquid holdup (unitless)
        
         %overall column performance - CO2
        Csat_CO2_in  = ppm0_CO2  * P / (1e6 * H_CO2);   %g/L at inlet gas conditions
        Csat_CO2_out = ppm_CO2_x(end) * P / (1e6 * H_CO2);  %g/L at outlet gas conditions
        kLa_L_CO2_all(f,j) = (1/max(t_x)) * log((Csat_CO2_in - c0_CO2)/(Csat_CO2_in - max(c_CO2_x))); %liquid-phase kLa 1/s
        kLa_T_CO2_all(f,j) = kLa_L_CO2_all(f,j) * liquid_holdup_tot; %total column kLa 1/s
        holdup_all(f,j) = liquid_holdup_tot;

        %overall column performance - O2
        Csat_O2_in  = ppm0_O2  * P / (1e6 * H_O2);   %g/L at inlet gas conditions
        Csat_O2_out = ppm_O2_x(end) * P / (1e6 * H_O2);  %g/L at outlet gas conditions
        kLa_L_O2_all(f,j) = (1/max(t_x)) * log((Csat_O2_in - c0_O2)/(Csat_O2_in - max(c_O2_x))); %liquid-phase kLa 1/s
        kLa_T_O2_all(f,j) = kLa_L_O2_all(f,j) * liquid_holdup_tot; %total column kLa 1/s
        holdup_all(f,j) = liquid_holdup_tot;
        
       %overall column performance - N2
        Csat_N2_in  = ppm0_N2  * P / (1e6 * H_N2);   %g/L at inlet gas conditions
        Csat_N2_out = ppm_N2_x(end) * P / (1e6 * H_N2);  %g/L at outlet gas conditions
        kLa_L_N2_all(f,j) = (1/max(t_x)) * log((Csat_N2_in - c0_N2)/(Csat_N2_in - max(c_N2_x))); %liquid-phase kLa 1/s
        kLa_T_N2_all(f,j) = kLa_L_N2_all(f,j) * liquid_holdup_tot; %total column kLa 1/s
        holdup_all(f,j) = liquid_holdup_tot;
        
        %CO2 purity
        total_dissolved = max(c_CO2_x)+max(c_O2_x)+max(c_N2_x);
        CO2_purity(f,j) = max(c_CO2_x)/total_dissolved*100;

        %capture efficiency
        if strcmp(flow_configs{f}, 'coflow')
            capture_eff_all(f,j) = (ppm0_CO2 - ppm_CO2_x(end)) / ppm0_CO2 * 100;
        else
            capture_eff_all(f,j) = (ppm0_CO2 - ppm_CO2_x(1))   / ppm0_CO2 * 100;
        end 
    end
end

%CO2 capture rate for each config and Ug
Q_g = Ug_range * A_c;  % m3/s
Csat_in      = ppm0_CO2 * P / (1e6 * H_CO2);   % g/L inlet equilibrium concentration
max_possible = Q * Csat_in * 1000;              % g/s liquid capacity limit

for f = 1:2
    CO2_mass_in          = Q_g * 1000 * Cg_t * ppm0_CO2/1e6;
    CO2_captured_all(f,:) = CO2_mass_in .* capture_eff_all(f,:)/100;
end

%--- Figure 1: Capture Rate and Efficiency ---
figure('color','white')

yyaxis left
h1 = semilogx(Ug_range, CO2_captured_all(1,:), 'o-',  'LineWidth', 2, 'Color', '#44AA99');
hold on
h2 = semilogx(Ug_range, CO2_captured_all(2,:), 'o--', 'LineWidth', 2, 'Color', '#44AA99');
h3 = yline(max_possible, ':', 'LineWidth', 1.5, 'Color', '#44AA99', ...
    'LabelHorizontalAlignment', 'left');
ylabel('CO_2 Capture Rate (g/s)')

yyaxis right
h4 = semilogx(Ug_range, capture_eff_all(1,:), 's-',  'LineWidth', 2, 'Color', '#CC6677');
h5 = semilogx(Ug_range, capture_eff_all(2,:), 's--', 'LineWidth', 2, 'Color', '#CC6677');
ylabel('CO_2 Capture Efficiency (%)')

xlabel('Gas Superficial Velocity, U_g (m/s)')
legend([h1, h2, h3, h4, h5], ...
    {'Capture Rate (co-flow)', 'Capture Rate (counter-flow)', ...
     'Liquid Saturation Limit', ...
     'Efficiency (co-flow)', 'Efficiency (counter-flow)'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 3)
fontsize(24,'points')

%--- Figure 2: CO2 Purity ---
figure('color','white')
semilogx(Ug_range, CO2_purity(1,:), 'o-',  'LineWidth', 2, 'Color', '#44AA99')
hold on
semilogx(Ug_range, CO2_purity(2,:), 'o--', 'LineWidth', 2, 'Color', '#44AA99')
xlabel('Gas Superficial Velocity, U_g (m/s)')
ylabel('CO_2 Purity in Dissolved Phase (%)')
legend({'Purity (co-flow)', 'Purity (counter-flow)'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal')
fontsize(24,'points')

%--- Figure 3: Dissolved concentration profiles at optimal Ug ---
% pick the shoulder point - approximately Ug = 0.1 m/s
[~, j_opt] = min(abs(Ug_range - 0.1));

figure('color','white')
subplot(1,2,1)
title('Co-flow')
hold on
plot(xspan, squeeze(c_CO2_out_all(1,j_opt,:)), 'LineWidth', 2, 'Color', '#44AA99', 'DisplayName', 'CO_2')
plot(xspan, squeeze(c_O2_out_all(1,j_opt,:)),  'LineWidth', 2, 'Color', '#CC6677', 'DisplayName', 'O_2')
plot(xspan, squeeze(c_N2_out_all(1,j_opt,:)),  'LineWidth', 2, 'Color', '#4477AA', 'DisplayName', 'N_2')
xlabel('Tower Height (m)')
ylabel('Dissolved Concentration (g/L)')
legend('Location','northwest')

subplot(1,2,2)
title('Counter-flow')
hold on
plot(xspan, squeeze(c_CO2_out_all(2,j_opt,:)), 'LineWidth', 2, 'Color', '#44AA99', 'DisplayName', 'CO_2')
plot(xspan, squeeze(c_O2_out_all(2,j_opt,:)),  'LineWidth', 2, 'Color', '#CC6677', 'DisplayName', 'O_2')
plot(xspan, squeeze(c_N2_out_all(2,j_opt,:)),  'LineWidth', 2, 'Color', '#4477AA', 'DisplayName', 'N_2')
xlabel('Tower Height (m)')
ylabel('Dissolved Concentration (g/L)')
legend('Location','northwest')
fontsize(24,'points')