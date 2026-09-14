%% 
clear all;close all;clc;

%DAC droplet tower model using Ionic Liquid [bmim][DCA]

%USER INPUTS---------------------------------------------------------------

%Import droplet size data OR Designate mean droplet diameter (mm) -
% select relavent input
DropletSize = readmatrix(['C:\Users\anbarron\Box\Research\Projects\DAC\Droplet Tower Modeling\DropletDiameterData.csv']); %droplet diameter data (mm)
%DropletSize = 1.6; %mean droplet diameter (mm)

%Input Tower Parameters
towerHeight = 15; %m
towerDiameter = 10; %m

%Input manifold parameters
d_orifice = 0.8/1000; %droplet manifold orifice diameter (m)
N_o = 30790750; %number of orifices in droplet manifold

%Input operating conditions
Q = 11.5465; %liquid flowrate (m3/s)
T = 21+273.15; %liquid temperature (K)

% Gas operating conditions - SET MASS FLOW FROM BASELINE
P_baseline = 20;  % atm - baseline design pressure (user input)
Ug_baseline = 0.025;  % m/s - desired Ug at baseline pressure (user input)

% For sequential/multistage - set operating pressure for THIS stage
P = 20; % atm - actual operating pressure (user input, = P_baseline for stage 1)
P_sat_sorbent = 0.0; %IL vapor pressure at setpoint T (atm)

%Relative Humidity 
RH = 0; % percent relatvie humidity 
Psat = 0.61078*exp(17.27*T/(T-273.15+237.3)); %saturation vapor pressure kPa
P_v = RH/100*Psat; %water vapor partial pressure kPa
P_d = P*101325 - P_v*1000; %air partial pressure Pa

%Gas Composition
MW_CO2 = 44.01; %CO2 gas molar mass (g/mol)
MW_O2 = 32; %O2 gas molar mass (g/mol)
MW_N2 = 28; %N2 gas molar mass (g/mol)
MW_Air = 28.97; %Air molar mass (g/mol)
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
R_d = 287.058; %ideal gas constant for dry air J/kg-K
R_v = 461.495; %ideal gas constant for water vapor J/kg-K
ppm0_CO2 = 400; %initial CO2 conentration (ppm)
ppm0_O2 = 210000; %initial O2 conentration (ppm)
ppm0_N2 = 780000; %initial N2 conentration (ppm)

c0_CO2 = 0; %initial CO2 dissolved gas concentration (g/L)
c0_O2 = 0; %initial O2 dissolved gas concentration (g/L)
c0_N2 = 0; %initial N2 dissolved gas concentration (g/L)

%Constant parameters
g = 9.8; %acceleration due to gravity (m/s2)
rho_g = P_d/(R_d*T)+(P_v*1000)/(R_v*T); %density of air (kg/m3)
mu_g = 1.86E-5; %dynamic viscocity of air (kg/m-s)

%IL Properties
MW_IL = 205.26; %g/mol [emim][DCA]
R_J = 8.314; %ideal gas constant (J/mol/K)
T_ref = 298.15; %K Reference temperatures for Arrhenius equations

%Temperature-Dependent IL Properties
rho_l = -0.55*T+1223; %IL density (kg/m3) DOI: 10.1021/je034261a
s_tension = (-0.0775E-3*T+71.88E-3)*1000; %surface tension (g/s2) DOI: 10.6100/IR639177
eta_l = 0.05417*exp(1000/(T-138)); %dynamic viscocity (cP) DOI: 10.6100/IR639177
DAB_CO2 = (6.7E5*(eta_l^-0.66)*(MW_IL^-0.89)*(rho_l/1000)^4.8*T^-3.3)/10000; %CO2 diffusivity (m2/s) DOI: 10.1021/ie070501u 
DAB_O2 = DAB_CO2; %approximation
DAB_N2 = DAB_CO2; %approximation

%Henry's Constant for CO2 DOI: 10.6100/IR639177
rho_l_CO2ref = -0.55*303+1223; %IL density (kg/m3) DOI: 10.1021/je034261a
Href_CO2 = 60.3*0.98*MW_IL/rho_l_CO2ref/MW_CO2; %reference H converted from bar to atm-L/g 
H_CO2 = Href_CO2*exp(-13020/R_J*(1/T-1/303)); %Henry's Constant (atm-L/g)

%Henry's Constant for O2 DOI: 10.1021/jp046404l - same as [bmim][Tf2N]
rho_l_O2ref = -0.55*323+1223; %IL density (kg/m3) DOI: 10.1021/je034261a
Href_O2 = 1110*0.98*MW_IL/rho_l_O2ref/MW_O2; %reference H converted from MPa to atm-L/g
H_O2 = Href_O2*exp(20600/R_J*(1/T-1/323)); %Henry's Constant (atm-L/g)

%Henry's Constant for N2 DOI: 10.1021/acs.iecr.5c03061
rho_l_N2ref = -0.55*308+1223; %IL density (kg/m3) DOI: 10.1021/je034261a
Href_N2 = 6000*0.98*MW_IL/rho_l_N2ref/MW_N2; %reference H converted from bar to atm-L/g
H_N2 = Href_N2; %Henry's Constant (atm-L/g) - enthalpy of absorption assumed 0 according to 10.1021/acs.iecr.5c03061 and 10.1021/acs.jced.3c00159

%END USER INPUTS ----------------------------------------------------------

% Derived geometry
A_orifice = pi*(d_orifice/2)^2 * N_o;   % m2
A_c = pi*(towerDiameter/2)^2;      % m2

% Compute mass flow from baseline conditions
Cg_t_baseline = P_baseline * MW_Air / (R * T); % g/m3 at baseline
Q_g_baseline = Ug_baseline * A_c;  % m3/s at baseline
mdot_gas = Q_g_baseline * Cg_t_baseline/1000; % kg/s - CONSERVED
mdot_liquid = Q*rho_l; %kg/s

% Compute actual operating conditions at stage pressure
Cg_t = P * MW_Air /(R*T);  % g/m3 at operating P
Q_g = mdot_gas / (Cg_t/1000); % m3/s at operating P
Ug_range = Q_g / A_c; % m/s at operating P
Lambda = 1e6 / (Ug_range * Cg_t/1000); % ppm*L/g/m - constant across stages


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
v0 = Q/A_orifice; %initial droplet velocity (m/s)
%time
t0 = 0; %inital time (s)

%intial conditions vector
y0 = [repmat([v0,c0_CO2,c0_O2,c0_N2,t0],1,M),ppm0_CO2,ppm0_O2,ppm0_N2];

%BVP solver options
options = bvpset('RelTol', 1e-4, 'AbsTol', 1e-6, 'NMax', 10000);

%Preallocate results - 3D: [flow config, Ug, x]
ppm_CO2_out_all  = zeros(1, length(Ug_range), length(xspan));
ppm_O2_out_all   = zeros(1, length(Ug_range), length(xspan));
ppm_N2_out_all   = zeros(1, length(Ug_range), length(xspan));
c_CO2_out_all    = zeros(1, length(Ug_range), length(xspan));
c_O2_out_all     = zeros(1, length(Ug_range), length(xspan));
c_N2_out_all     = zeros(1, length(Ug_range), length(xspan));
v_out_all        = zeros(1, length(Ug_range), length(xspan));
t_out_all        = zeros(1, length(Ug_range), length(xspan));
capture_eff_all  = zeros(1, length(Ug_range));
CO2_purity       = zeros(1, length(Ug_range));
kLa_L_CO2_all   = zeros(1, length(Ug_range));
kLa_T_CO2_all   = zeros(1, length(Ug_range));
kLa_L_O2_all    = zeros(1, length(Ug_range));
kLa_T_O2_all    = zeros(1, length(Ug_range));
kLa_L_N2_all    = zeros(1, length(Ug_range));
kLa_T_N2_all    = zeros(1, length(Ug_range));
holdup_all       = zeros(1, length(Ug_range));
CO2_captured_all = zeros(1, length(Ug_range));


for j = 1:length(Ug_range)

    Ug = Ug_range(j);
    lambda = 10^6/(Ug*Cg_t/1000); %conversion from mdot to ppm

    constants = [g, rho_g, rho_l, mu_g, DAB_CO2, DAB_O2, DAB_N2, s_tension, P, H_CO2, H_O2, H_N2, A_c, Q, lambda];

    %counter-flow - BVP solved with bvp4c
    if j == 1
        solinit = bvpinit(xspan, y0);
    else
        solinit = bvpinit(xspan, @(x) interp1(xspan, y_out, x)');
    end
    sol   = bvp4c(@(x,y) ODEsys_Angelo_CounterFlowcoupled(x, y, GQMOM_model_x, GQMOM_model_w, constants), @(ya,yb) BC_coupled(ya, yb, v0, c0_CO2,c0_O2,c0_N2, ppm0_CO2, ppm0_O2, ppm0_N2, M),solinit, options);
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
    ppm_CO2_out_all(1,j,:) = ppm_CO2_x;
    ppm_O2_out_all(1,j,:) = ppm_O2_x;
    ppm_N2_out_all(1,j,:) = ppm_N2_x;
    c_CO2_out_all(1,j,:)   = c_CO2_x;
    c_O2_out_all(1,j,:)   = c_O2_x;
    c_N2_out_all(1,j,:)   = c_N2_x;
    v_out_all(1,j,:)   = v_x;
    t_out_all(1,j,:)   = t_x;

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
    Csat_CO2_out = ppm_CO2_x(1) * P / (1e6 * H_CO2);  %g/L at outlet gas conditions
    kLa_L_CO2_all(1,j) = (1/max(t_x)) * log((Csat_CO2_in - c0_CO2)/(Csat_CO2_in - max(c_CO2_x))); %liquid-phase kLa 1/s
    kLa_T_CO2_all(1,j) = kLa_L_CO2_all(1,j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(1,j) = liquid_holdup_tot;

    %overall column performance - O2
    Csat_O2_in  = ppm0_O2  * P / (1e6 * H_O2);   %g/L at inlet gas conditions
    Csat_O2_out = ppm_O2_x(1) * P / (1e6 * H_O2);  %g/L at outlet gas conditions
    kLa_L_O2_all(1,j) = (1/max(t_x)) * log((Csat_O2_in - c0_O2)/(Csat_O2_in - max(c_O2_x))); %liquid-phase kLa 1/s
    kLa_T_O2_all(1,j) = kLa_L_O2_all(1,j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(1,j) = liquid_holdup_tot;

    %overall column performance - N2
    Csat_N2_in  = ppm0_N2  * P / (1e6 * H_N2);   %g/L at inlet gas conditions
    Csat_N2_out = ppm_N2_x(1) * P / (1e6 * H_N2);  %g/L at outlet gas conditions
    kLa_L_N2_all(1,j) = (1/max(t_x)) * log((Csat_N2_in - c0_N2)/(Csat_N2_in - max(c_N2_x))); %liquid-phase kLa 1/s
    kLa_T_N2_all(1,j) = kLa_L_N2_all(1,j) * liquid_holdup_tot; %total column kLa 1/s
    holdup_all(1,j) = liquid_holdup_tot;

    %CO2 purity
    total_dissolved = max(c_CO2_x)+max(c_O2_x)+max(c_N2_x);
    CO2_purity(1,j) = max(c_CO2_x)/total_dissolved*100;

    %capture efficiency
    capture_eff_CO2(1,j) = (ppm0_CO2 - ppm_CO2_x(1))   / ppm0_CO2 * 100;
    capture_eff_O2(1,j) = (ppm0_O2 - ppm_O2_x(1))   / ppm0_O2 * 100;
    capture_eff_N2(1,j) = (ppm0_N2 - ppm_N2_x(1))   / ppm0_N2 * 100;

end


%CO2 capture rate for each config and Ug
Q_g = Ug_range * A_c;  % m3/s
Csat_in      = ppm0_CO2 * P / (1e6 * H_CO2);   % g/L inlet equilibrium concentration
max_possible = Q * Csat_in * 1000;              % g/s liquid capacity limit


for f = 1
    CO2_mass_in = Q_g*Cg_t*ppm0_CO2/1e6;
    CO2_captured_all(f,:) = CO2_mass_in .* capture_eff_CO2(f,:)/100;
end

%gas capture rate 
CO2_capture_rate = (max(c_CO2_x)-min(c_CO2_x))*Q*1000; %g/s
O2_capture_rate = (max(c_O2_x)-min(c_O2_x))*Q*1000; %g/s
N2_capture_rate = (max(c_N2_x)-min(c_N2_x))*Q*1000; %g/s

%--- Figure 1: Capture Rate and Efficiency ---
figure('color','white')
ax=gca;

yyaxis left
h1 = semilogx(Ug_range, CO2_captured_all(1,:), 'o-',  'LineWidth', 2, 'Color', '#44AA99');
hold on
h3 = yline(max_possible, ':', 'LineWidth', 1.5, 'Color', '#44AA99', ...
    'LabelHorizontalAlignment', 'left');
ylabel('CO_2 Capture Rate (g/s)','Color', '#44AA99');
ax.YAxis(1).Color='#44AA99';

yyaxis right
h4 = semilogx(Ug_range, capture_eff_CO2(1,:), 's-',  'LineWidth', 2, 'Color', '#CC6677');
ylabel('CO_2 Capture Efficiency (%)','Color', '#CC6677');
ax.YAxis(2).Color='#CC6677';

xlabel('Gas Superficial Velocity, U_g (m/s)')
legend([h1, h3, h4], ...
    {'Capture Rate', ...
    'Liquid Saturation Limit', ...
    'Efficiency'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 3)
fontsize(24,'points')

%--- Figure 2: CO2 Purity ---
figure('color','white')
semilogx(Ug_range, CO2_purity(1,:), 'o-',  'LineWidth', 2, 'Color', '#44AA99')
hold on
xlabel('Gas Superficial Velocity, U_g (m/s)')
ylabel('CO_2 Purity in Dissolved Phase (%)')
fontsize(24,'points')

%--- Figure 3: Dissolved concentration profiles at optimal Ug ---
% pick the shoulder point - approximately Ug = 0.05 m/s
[~, j_opt] = min(abs(Ug_range - 0.05));

figure('color','white')
hold on
plot(xspan, squeeze(c_CO2_out_all(1,j_opt,:)), 'LineWidth', 2, 'Color', '#44AA99', 'DisplayName', 'CO_2')
plot(xspan, squeeze(c_O2_out_all(1,j_opt,:)),  'LineWidth', 2, 'Color', '#CC6677', 'DisplayName', 'O_2')
plot(xspan, squeeze(c_N2_out_all(1,j_opt,:)),  'LineWidth', 2, 'Color', '#4477AA', 'DisplayName', 'N_2')
xlabel('Tower Height (m)')
ylabel('Dissolved Concentration (g/L)')
legend('Location','northwest')
fontsize(24,'points')

%exporting results 
results.IL_name     = '[bmim][DCA]';
results.Ug_range    = Ug_range;

% absorption results
results.CO2_captured    = CO2_captured_all;
results.capture_eff     = capture_eff_all;
results.CO2_purity_abs  = CO2_purity;

% profiles at optimal Ug
[~, j_opt] = min(abs(Ug_range - 0.05));
results.c_CO2_profile   = squeeze(c_CO2_out_all(1,j_opt,:));  % counter-flow
results.c_O2_profile    = squeeze(c_O2_out_all(1,j_opt,:));
results.c_N2_profile    = squeeze(c_N2_out_all(1,j_opt,:));
results.ppm_CO2_profile = squeeze(ppm_CO2_out_all(1,j_opt,:));

