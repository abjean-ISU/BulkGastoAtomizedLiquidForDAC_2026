clc; clear all; close all;

%% User Inputs

%Droplet Data
DropletSize = readmatrix(['C:\Users\anbarron\Box\Research\Projects\DAC\Droplet Tower Modeling\DropletDiameterData.csv']); %droplet diameter data (mm)

%Input Tower Parameters
towerDiameter = 7; %m
A_c = pi*(towerDiameter/2)^2; % m2

%Input manifold parameters
d_orifice = 0.8/1000; %droplet manifold orifice diameter (m)
N_o = 373928; %number of orifices in droplet manifold

%Input operating conditions
Q = 0.1402; %liquid flowrate (m3/s)
T = 25+273.15; %liquid temperature (K)
P_set = 3; %Operating pressure setpoint (atm) 

%Relative Humidity 
RH = 0; % percent relatvie humidity 
Psat = 0.61078*exp(17.27*T/(T-273.15+237.3)); %saturation vapor pressure kPa
P_v = RH/100*Psat; %water vapor partial pressure kPa
P_d = P_set*101325 - P_v*1000; %air partial pressure Pa

%Constants 
g = 9.8; %acceleration due to gravity (m/s2)
R_d = 287.058; %ideal gas constant for dry air J/kg-K
R_v = 461.495; %ideal gas constant for water vapor J/kg-K
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
MW_Air = 28.97; %Air molar mass (g/mol)

% Gas baseline operating conditions - SET MASS FLOW FROM BASELINE
P_baseline = 10;  % atm - baseline design pressure (user input)
Ug_baseline = 0.05;  % m/s - desired Ug at baseline pressure (user input)
rho_g_baseline = P_baseline*101325/(R_d*T); %baseline gas sensity 
mdot_gas = rho_g_baseline * Ug_baseline * A_c; % kg/s

%Gas Properties at setpoint 
mu_g = 1.86E-5; %dynamic viscocity of air (kg/m-s)
rho_g = P_d/(R_d*T)+(P_v*1000)/(R_v*T); %density of air (kg/m3)
Cg_t = P_set * MW_Air /(R*T);  % g/m3 at operating P
Q_g = mdot_gas / (Cg_t/1000); % m3/s at operating P
Ug = Q_g / A_c; % m/s at operating P

%liquid properties 
rho_l = 997; %density of water (kg/m3)

%% end user inputs 


r_drops = (DropletSize/1000)/2; % mm diameter -> m radius
vol_drops = (4/3)*pi*r_drops.^ 3; % volume per droplet
vol_total = sum(vol_drops);

M = length(r_drops); %number of droplets

%-- Terminal velocity for each quadrature point -----------------------
vt = zeros(M, 1);

for i = 1:M
    r = r_drops(i);

    % Stokes initial guess
    vt_guess = (2/9) * (rho_l - rho_g) * g * r^2 / mu_g;
    vt_guess = max(vt_guess, 1e-6);

    % Fixed-point iteration
    for iter = 1:200
        Re_t  = 2 * rho_g * vt_guess * r / mu_g;  
        accel = g - (9/2) * (mu_g / (r^2 * rho_l)) * vt_guess * (1 + 0.158 * Re_t^(2/3));
        vt_new = vt_guess + accel * 0.005;
        if abs(vt_new - vt_guess) < 1e-10
            break
        end
        vt_guess = max(vt_new, 1e-6);
    end

    vt(i) = vt_guess;
end

%-- Entrainment criterion ---------------------------------------------
entrained = vt < Ug;
vol_frac_entrained = sum(vol_drops(entrained))/vol_total*100;


%-- Critical diameter -------------------------------------------------
if any(entrained)
    d_critical_mm = 2 * max(r_drops(entrained))*1000;
else
    d_critical_mm = 0;
end

