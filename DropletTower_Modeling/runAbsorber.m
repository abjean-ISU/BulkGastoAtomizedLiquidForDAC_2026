function [out, y_out] = runAbsorber(c, p, GQMOM_model_w, GQMOM_model_x, y_guess)
%counter-flow droplet absorber for one case
%c: case row from BGAL_cases.xlsx; p: sorbent properties (props_water / props_IL)
%y_guess: previous solution on the same xspan for the bvp4c initial guess ([] for none)

%constants
g = 9.8; %acceleration due to gravity (m/s2)
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
R_d = 287.058; %ideal gas constant for dry air J/kg-K
MW_CO2 = 44.01; %CO2 gas molar mass (g/mol)
MW_O2 = 32; %O2 gas molar mass (g/mol)
MW_N2 = 28; %N2 gas molar mass (g/mol)
MW_Air = 28.97; %Air molar mass (g/mol)
mu_g = 1.86E-5; %dynamic viscocity of air (kg/m-s)

%case inputs
T = c.T_C+273.15; %liquid temperature (K)
P = c.P; %operating pressure (atm)
Q = c.Q; %liquid flowrate (m3/s)
A_c = c.A_c; %column cross sectional area (m2)
A_orifice = pi*(c.d_orifice_mm/1000/2)^2*c.N_o; %total orifice area (m2)
rho_g = P*101325/(R_d*T); %density of dry air (kg/m3)
H_CO2 = p.H_CO2*c.H_scale; %Henry's Constants (atm-L/g)
H_O2 = p.H_O2*c.H_scale;
H_N2 = p.H_N2*c.H_scale;

%gas mass flow set at baseline conditions, Ug at operating pressure
Cg_t_baseline = c.P_baseline*MW_Air/(R*T); %g/m3 at baseline
mdot_gas = c.Ug_baseline*A_c*Cg_t_baseline/1000; %kg/s
Cg_t = P*MW_Air/(R*T); %g/m3 at operating P
Q_g = mdot_gas/(Cg_t/1000); %m3/s
Ug = Q_g/A_c; %m/s
lambda = 1e6/(Ug*P/(R*T)); %molar gas balance, same for all gases

M = length(GQMOM_model_x); %number of quadrature nodes
xspan = 0:0.01:c.towerHeight; %location (m)
v0 = Q/A_orifice; %initial droplet velocity (m/s)
y0 = [repmat([v0,c.c0_CO2,c.c0_O2,c.c0_N2,0],1,M),c.ppm0_CO2,c.ppm0_O2,c.ppm0_N2];

constants = [g, rho_g, p.rho_l, mu_g, p.DAB_CO2, p.DAB_O2, p.DAB_N2, p.s_tension, P, H_CO2, H_O2, H_N2, A_c, Q, lambda, MW_CO2, MW_O2, MW_N2];

%counter-flow - BVP solved with bvp4c
options = bvpset('RelTol', 1e-4, 'AbsTol', 1e-6, 'NMax', 10000);
if ~isempty(y_guess) && size(y_guess,1) == length(xspan)
    solinit = bvpinit(xspan, @(x) interp1(xspan, y_guess, x)');
else
    solinit = bvpinit(xspan, y0);
end
sol = bvp4c(@(x,y) ODEsys_absorber(x, y, GQMOM_model_x, GQMOM_model_w, constants), ...
    @(ya,yb) BC_coupled(ya, yb, v0, c.c0_CO2, c.c0_O2, c.c0_N2, c.ppm0_CO2, c.ppm0_O2, c.ppm0_N2, M), solinit, options);
y_out = deval(sol, xspan)';

for i = 1:M
    idx = (i-1)*5 + 1;
    v_ri(:,i) = y_out(:,idx);
    c_CO2_ri(:,i) = y_out(:,idx+1);
    c_O2_ri(:,i) = y_out(:,idx+2);
    c_N2_ri(:,i) = y_out(:,idx+3);
    t_ri(:,i) = y_out(:,idx+4);
end

%QMOM averages: velocity number-weighted, concentration and time volume-weighted
volW = GQMOM_model_w(:).*GQMOM_model_x(:).^3; volW = volW/sum(volW);
v_x = v_ri*GQMOM_model_w(:);
c_CO2_x = c_CO2_ri*volW;
c_O2_x = c_O2_ri*volW;
c_N2_x = c_N2_ri*volW;
t_x = t_ri*volW;
ppm_CO2_x = y_out(:,end-2);
ppm_O2_x = y_out(:,end-1);
ppm_N2_x = y_out(:,end);

%liquid holdup
tower_v_tot = A_c*c.towerHeight; %tower volume (m3)
droplet_v = sum(GQMOM_model_w(:).*(4/3)*pi.*GQMOM_model_x(:).^3); %mean droplet volume (m3)
droplet_v_tot_m3 = Q*t_x(end); %total droplet (LH) volume (m3)
N_tot = droplet_v_tot_m3/droplet_v; %total number of droplets
liquid_holdup_tot = droplet_v_tot_m3/tower_v_tot; %total liquid holdup (unitless)

%liquid-phase kLa (1/s)
Csat_CO2_in = c.ppm0_CO2*P/(1e6*H_CO2); %g/L at inlet gas conditions
Csat_O2_in = c.ppm0_O2*P/(1e6*H_O2);
Csat_N2_in = c.ppm0_N2*P/(1e6*H_N2);
kLa_L_CO2 = (1/t_x(end))*log((Csat_CO2_in-c.c0_CO2)/(Csat_CO2_in-c_CO2_x(end)));
kLa_L_O2 = (1/t_x(end))*log((Csat_O2_in-c.c0_O2)/(Csat_O2_in-c_O2_x(end)));
kLa_L_N2 = (1/t_x(end))*log((Csat_N2_in-c.c0_N2)/(Csat_N2_in-c_N2_x(end)));

%capture efficiency (gas side, mole ppm)
capture_eff_CO2 = (c.ppm0_CO2-ppm_CO2_x(1))/c.ppm0_CO2*100;
capture_eff_O2 = (c.ppm0_O2-ppm_O2_x(1))/c.ppm0_O2*100;
capture_eff_N2 = (c.ppm0_N2-ppm_N2_x(1))/c.ppm0_N2*100;

%CO2 capture rate (g/s): gas side vs liquid side, and liquid capacity limit
CO2_mass_in = Q_g*P*MW_CO2/(R*T)*c.ppm0_CO2/1e6;
CO2_captured = CO2_mass_in*capture_eff_CO2/100;
CO2_capture_rate = (c_CO2_x(end)-c_CO2_x(1))*Q*1000;
max_possible = Q*Csat_CO2_in*1000;

%CO2 purity in dissolved gas (mass basis)
CO2_purity = c_CO2_x(end)/(c_CO2_x(end)+c_O2_x(end)+c_N2_x(end))*100;

out = struct('caseName', c.caseName, 'stage', 'abs', 'sorbent', c.sorbent, 'Ug', Ug, 'mdot_gas', mdot_gas, ...
    'mdot_liquid', Q*p.rho_l, 'eff_CO2', capture_eff_CO2, 'eff_O2', capture_eff_O2, 'eff_N2', capture_eff_N2, ...
    'CO2_captured', CO2_captured, 'CO2_capture_rate', CO2_capture_rate, 'max_possible', max_possible, ...
    'CO2_purity', CO2_purity, 'liquid_holdup_tot', liquid_holdup_tot, 'N_tot', N_tot, ...
    'kLa_L_CO2', kLa_L_CO2, 'kLa_L_O2', kLa_L_O2, 'kLa_L_N2', kLa_L_N2, 'kLa_T_CO2', kLa_L_CO2*liquid_holdup_tot, ...
    'xspan', xspan(:), 'v_x', v_x, 't_x', t_x, 'c_CO2_x', c_CO2_x, 'c_O2_x', c_O2_x, 'c_N2_x', c_N2_x, ...
    'ppm_CO2_x', ppm_CO2_x, 'ppm_O2_x', ppm_O2_x, 'ppm_N2_x', ppm_N2_x, 'props', p, 'inputs', c);

end
