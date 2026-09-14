clear all;close all;clc;

%DAC droplet tower model 

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
DAB = 1.67E-9; %diffusion coefficient of CO2 in water at correct Temp (m2/s)
c0 = 6.35E-04; %initial dissolved gas concentration (g/L)
mass_Gas = 44.01; %CO2 gas molar mass (g/mol)
mass_Air = 28.97; %Air molar mass (g/mol)
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
ppm0 = 3.80E+02; %initial CO2 conentration (ppm)
Ug_range = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0];  % m/s



%END USER INPUTS ----------------------------------------------------------

Cg_t = P*mass_Air/(R*T)/1000; %total concentraiton of air (g/L)

H = (1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T)))/55324*1000/mass_Gas; %Henry's Constant for CO2 (atm-L/g)

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
y0 = [repmat([v0,c0,t0],1,M),ppm0];

% Preallocate results storage
ppm_out        = zeros(length(Ug_range), length(xspan));
c_out          = zeros(length(Ug_range), length(xspan));
v_out          = zeros(length(Ug_range), length(xspan));
t_out          = zeros(length(Ug_range), length(xspan));
capture_eff    = zeros(length(Ug_range), 1);
kLa_L_all      = zeros(length(Ug_range), 1);
kLa_T_all      = zeros(length(Ug_range), 1);
holdup_all     = zeros(length(Ug_range), 1);

for j = 1:length(Ug_range)
    
    Ug = Ug_range(j); 
    lambda = 10^6/(Ug*Cg_t); %conversion from mdot to ppm

    %Constant parameters 
    g = 9.8; %acceleration due to gravity (m/s2)
    rho_g = 1.184; %density of air (kg/m3)
    rho_l = 997; %density of water (kg/m3) 
    mu_g = 1.53E-5; %viscocity of air (m2/s)
    s_tension = 72.8; %water surface tenstion (g/s2)
    
    constants = [g, rho_g, rho_l, mu_g, DAB, s_tension, P, H, A_c, Q, lambda]; 
    
    %Gas absorption modeling <v(x|r)> and <c(x|r)> using kL model by Angelo et al.
    [x_out,y_out] = ode45(@(x,y) ODEsys_Angelo_CoFlowcoupled(x,y,GQMOM_model_x,GQMOM_model_w,constants), xspan, y0);
    
    %extracting size class results
    for i=1:M
        idx = (i-1)*3 + 1;
        v_ri(:,i) = y_out(:,idx);
        c_ri(:,i) = y_out(:,idx+1);
        t_ri(:,i) = y_out(:,idx+2);
    end
    
    
    %QMOM weighted averages
    v_x = sum(v_ri .* GQMOM_model_w', 2);
    c_x = sum(c_ri .* GQMOM_model_w', 2);
    t_x = sum(t_ri .* GQMOM_model_w', 2);
    ppm_x = y_out(:,end);  %gas phase CO2 concentration (ppm)

    %store profiles 
    ppm_out(j,:) = ppm_x;
    c_out(j,:) = c_x;
    v_out(j,:) = v_x;
    t_out(j,:) = t_x;
    
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
    
    %overall column performance 
    Csat_in  = ppm0  * P / (1e6 * H);   %g/L at inlet gas conditions
    Csat_out = ppm_x(end) * P / (1e6 * H);  %g/L at outlet gas conditions
    kLa_L_all(j) = (1/max(t_x)) * log((Csat_in - c0)/(Csat_in - max(c_x))); %liquid-phase kLa 1/s
    kLa_T_all(j) = kLa_L_all(j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(j) = liquid_holdup_tot;
    
    %capture efficienty
    capture_eff(j) = (ppm0 - ppm_x(end))/ppm0*100; %percent change 
end

% total CO2 capture rate for each Ug
Q_g = Ug_range * A_c; % m3/s gas flowrate
CO2_mass_in = Q_g *1000 * Cg_t * ppm0/1e6; % g/s CO2 entering
CO2_captured = CO2_mass_in .* capture_eff'/100; % g/s CO2 captured

% maximum possible capture if droplets reach equilibrium
Csat_in = ppm0 * P / (1e6 * H); % g/L equilibrium conc
max_possible = Q * Csat_in * 1000; % g/s - liquid capacity limit

figure('color','white')

% left axis - capture rate
yyaxis left
semilogx(Ug_range, CO2_captured, 'o-', 'LineWidth', 2, 'Color', '#44AA99')
ylabel('CO_2 Capture Rate (g/s)')
yline(max_possible, '--', 'LineWidth', 1.5, 'Color', '#44AA99', ...
    'Label', 'Liquid saturation limit', 'LabelHorizontalAlignment', 'left')

% right axis - capture efficiency
yyaxis right
semilogx(Ug_range, capture_eff, 's--', 'LineWidth', 2, 'Color', '#CC6677')
ylabel('CO_2 Capture Efficiency (%)')

% formatting
xlabel('Gas Superficial Velocity, U_g (m/s)')
legend({'Capture Rate', 'Liquid Saturation Limit', 'Capture Efficiency'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal')
fontsize(24, 'points')