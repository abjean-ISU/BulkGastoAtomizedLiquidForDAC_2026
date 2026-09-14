clear all;close all;clc;

%DAC droplet tower model with water

%USER INPUTS---------------------------------------------------------------

%Import droplet size data OR Designate mean droplet diameter (mm) -
% select relavent input
DropletSize = readmatrix(['C:\Users\anbarron\Box\Research\Projects\DAC\Droplet Tower Modeling\DropletDiameterData.csv']); %droplet diameter data (mm)
%DropletSize = 1.6; %mean droplet diameter (mm)

%Input Tower Parameters
towerHeight = 25; %m
towerDiameter = 10; %m

%Input manifold parameters
d_orifice = 0.8/1000; %droplet manifold orifice diameter (m)
N_o = 1539537; %number of orifices in droplet manifold

%Input operating conditions
Q = 0.5773; %liquid flowrate (m3/s)
T = 25+273.15; %liquid temperature (K)
P_sat_sorbent = 0.03; %water vapor pressure at setpoint T (atm)
P = 0.05;  %operating (headspace) pressure (atm)

%Gas Composition
DAB_CO2 = 1.67E-9; %diffusion coefficient of CO2 in water at correct Temp (m2/s)
DAB_O2 = 2.01E-9; %diffusion coefficient of O2
DAB_N2 = 2.00E-9; %diffusion coefficient of N2
MW_CO2 = 44.01; %CO2 gas molar mass (g/mol)
MW_O2 = 32; %O2 gas molar mass (g/mol)
MW_N2 = 28; %N2 gas molar mass (g/mol)
MW_Air = 28.97; %Air molar mass (g/mol)
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
ppm0_CO2 = 0.0001; %initial CO2 conentration (ppm)
ppm0_O2 = 0.0001; %initial O2 conentration (ppm)
ppm0_N2 = 0.0001; %initial N2 conentration (ppm)
ppm_vac = 0.0001;

c0_CO2 =0.0166; %initial CO2 dissolved gas concentration (g/L)
c0_O2 = 0.0015; %initial O2 dissolved gas concentration (g/L)
c0_N2 = 0.0028; %initial N2 dissolved gas concentration (g/L)

%Constant parameters
g = 9.8; %acceleration due to gravity (m/s2)
rho_l = 997; %density of water (kg/m3)
mu_g = 1.86E-5; %dynamic viscocity of air (kg/m-s)
s_tension = 72.8; %water surface tenstion (g/s2)

%END USER INPUTS ----------------------------------------------------------

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
options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NMax', 10000);

%Preallocate results - 3D: [flow config, Ug, x]
ppm_CO2_out_all  = zeros(1, length(P), length(xspan));
ppm_O2_out_all   = zeros(1, length(P), length(xspan));
ppm_N2_out_all   = zeros(1, length(P), length(xspan));
c_CO2_out_all    = zeros(1, length(P), length(xspan));
c_O2_out_all     = zeros(1, length(P), length(xspan));
c_N2_out_all     = zeros(1, length(P), length(xspan));
v_out_all        = zeros(1, length(P), length(xspan));
t_out_all        = zeros(1, length(P), length(xspan));
desorption_eff_all  = zeros(1, length(P));
CO2_product_purity       = zeros(1, length(P));
kLa_L_CO2_all   = zeros(1, length(P));
kLa_T_CO2_all   = zeros(1, length(P));
kLa_L_O2_all    = zeros(1, length(P));
kLa_T_O2_all    = zeros(1, length(P));
kLa_L_N2_all    = zeros(1, length(P));
kLa_T_N2_all    = zeros(1, length(P));
holdup_all       = zeros(1, length(P));

for j = 1:length(P)

    Cg_t = P*MW_Air/(R*T)/1000; %total concentraiton of air (g/L)
    rho_g = P*MW_Air/(R*T)/1000; 

    constants = [g, rho_g, rho_l, mu_g, DAB_CO2, DAB_O2, DAB_N2, s_tension, P, H_CO2, H_O2, H_N2, A_c, Q, ppm_vac];


    %counter-flow - BVP solved with bvp4c
    if j == 1
        solinit = bvpinit(xspan, y0);
    else
        solinit = bvpinit(xspan, @(x) interp1(xspan, y_out, x)');
    end
    sol   = bvp4c(@(x,y) ODEsys_Angelo_Desorption(x, y, GQMOM_model_x, GQMOM_model_w, constants), @(ya,yb) BC_coupled_desorption(ya, yb, v0, c0_CO2, c0_O2, c0_N2, ppm_vac, M),solinit, options);
    y_out = deval(sol, xspan)';

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
    ppm_CO2_out_all(1,j,:) = ppm_CO2_x';
    ppm_O2_out_all(1,j,:) = ppm_O2_x';
    ppm_N2_out_all(1,j,:) = ppm_N2_x';
    c_CO2_out_all(1,j,:)   = c_CO2_x';
    c_O2_out_all(1,j,:)   = c_O2_x';
    c_N2_out_all(1,j,:)   = c_N2_x';
    v_out_all(1,j,:)   = v_x';
    t_out_all(1,j,:)   = t_x';

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
    kLa_L_CO2_all(1,j) = (1/max(t_x)) * log((c0_CO2 - Csat_CO2_in)/(min(c_CO2_x)-Csat_CO2_in)); %liquid-phase kLa 1/s
    kLa_T_CO2_all(1,j) = kLa_L_CO2_all(1,j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(1,j) = liquid_holdup_tot;

    %overall column performance - O2
    Csat_O2_in  = ppm0_O2  * P / (1e6 * H_O2);   %g/L at inlet gas conditions
    Csat_O2_out = ppm_O2_x(end) * P / (1e6 * H_O2);  %g/L at outlet gas conditions
    kLa_L_O2_all(1,j) = (1/max(t_x)) * log((c0_O2 - Csat_O2_in)/(min(c_O2_x)-Csat_O2_in)); %liquid-phase kLa 1/s
    kLa_T_O2_all(1,j) = kLa_L_O2_all(1,j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(1,j) = liquid_holdup_tot;

    %overall column performance - N2
    Csat_N2_in  = ppm0_N2  * P / (1e6 * H_N2);   %g/L at inlet gas conditions
    Csat_N2_out = ppm_N2_x(end) * P / (1e6 * H_N2);  %g/L at outlet gas conditions
    kLa_L_N2_all(1,j) = (1/max(t_x)) * log((c0_N2 - Csat_N2_in)/(min(c_N2_x)-Csat_N2_in)); %liquid-phase kLa 1/s
    kLa_T_N2_all(1,j) = kLa_L_N2_all(1,j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(1,j) = liquid_holdup_tot;

    %CO2 purity
    ppm_CO2_out = ppm_CO2_x(1);
    ppm_O2_out = ppm_O2_x(1);
    ppm_N2_out = ppm_N2_x(1);
    CO2_product_purity(1,j) = ppm_CO2_out/(ppm_CO2_out+ppm_O2_out+ppm_N2_out)*100;

    %capture efficiency
    desorption_eff_CO2(1,j) = (c0_CO2 - c_CO2_x(end)) / c0_CO2 * 100;
    desorption_eff_O2(1,j) = (c0_O2 - c_O2_x(end)) / c0_O2 * 100;
    desorption_eff_N2(1,j) = (c0_N2 - c_N2_x(end)) / c0_N2 * 100;

end


%CO2 desorption rate for each config and P
CO2_mass_released(1,:) = Q * c0_CO2 .* desorption_eff_CO2(1,:)/100 * 1000;  % g/s

%gas capture rate 
CO2_desorption_rate = (c_CO2_x(1)-c_CO2_x(end))*Q*1000; %g/s
O2_desorption_rate = (c_O2_x(1)-c_O2_x(end))*Q*1000; %g/s
N2_desorption_rate = (c_N2_x(1)-c_N2_x(end))*Q*1000; %g/s

CO2_des_mol = CO2_desorption_rate / MW_CO2;  % mol/s
O2_des_mol  = O2_desorption_rate  / MW_O2;   % mol/s
N2_des_mol  = N2_desorption_rate  / MW_N2;   % mol/s

% Water vapor from desorption efficiency
E_CO2_actual = (c_CO2_x(1) - c_CO2_x(end)) / c_CO2_x(1);
y_H2O_actual = E_CO2_actual * P_sat_sorbent / P;
if y_H2O_actual >= 1, y_H2O_actual = 0.99; end

F_gas_mol = CO2_des_mol + O2_des_mol + N2_des_mol;
H2O_evap_mol = F_gas_mol * y_H2O_actual / (1 - y_H2O_actual);  % mol/s

% Actual Q_g at vacuum conditions
Q_g = (F_gas_mol + H2O_evap_mol) * R * T / P;  % m3/s
gas_mas = (CO2_desorption_rate+O2_desorption_rate+N2_desorption_rate+H2O_evap_mol*18)/1000; %kg/s
Ug_actual = Q_g / A_c;  % m/s — for reporting

%--- Figure 1: Desorption Rate and Efficiency ---
figure('color','white')
ax=gca;

yyaxis left
h1 = semilogx(P, CO2_mass_released(1,:), 'o-',  'LineWidth', 2, 'Color', '#44AA99');
hold on
ylabel('CO_2 Desorption Rate (g/s)')
ax.YAxis(1).Color='#44AA99';

yyaxis right
h3 = semilogx(P, desorption_eff_CO2(1,:), 's-',  'LineWidth', 2, 'Color', '#CC6677');
hold on
ylabel('CO_2 Desorption Efficiency (%)')
ax.YAxis(2).Color='#CC6677';

xlabel('Vacuum Pressure, P (atm)')
legend([h1, h3], ...
    { 'Desorption Rate', ...
    'Efficiency'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 2)
fontsize(24,'points')

%--- Figure 2: Dissolved concentration profiles at optimal P ---
% pick lowest pressure as most aggressive desorption condition
[~, j_opt] = min(P);

figure('color','white')
hold on
plot(xspan, squeeze(c_CO2_out_all(1,j_opt,:)), 'LineWidth', 2, 'Color', '#44AA99', 'DisplayName', 'CO_2')
plot(xspan, squeeze(c_O2_out_all(1,j_opt,:)),  'LineWidth', 2, 'Color', '#CC6677', 'DisplayName', 'O_2')
plot(xspan, squeeze(c_N2_out_all(1,j_opt,:)),  'LineWidth', 2, 'Color', '#4477AA', 'DisplayName', 'N_2')
yline(c0_CO2, '--', 'Color', '#44AA99', 'LineWidth', 1, 'DisplayName', 'C_{0,CO2}')
yline(c0_O2,  '--', 'Color', '#CC6677', 'LineWidth', 1, 'DisplayName', 'C_{0,O2}')
yline(c0_N2,  '--', 'Color', '#4477AA', 'LineWidth', 1, 'DisplayName', 'C_{0,N2}')
xlabel('Tower Height (m)')
ylabel('Dissolved Concentration (g/L)')
legend('Location','northeast')
fontsize(24,'points')
