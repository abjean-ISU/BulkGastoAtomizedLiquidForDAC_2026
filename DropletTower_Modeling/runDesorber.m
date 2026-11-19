function [out, y_out] = runDesorber(c, p, GQMOM_model_w, GQMOM_model_x, y_guess)
%vacuum droplet desorber for one case
%c: case row from BGAL_cases.xlsx; p: sorbent properties (props_water / props_IL)
%y_guess: previous solution on the same xspan for the bvp4c initial guess ([] for none)

%constants
g = 9.8; %acceleration due to gravity (m/s2)
R = 8.205e-5; %ideal gas constant m3-atm/mol-K
MW_Air = 28.97; %Air molar mass (g/mol)
mu_g = 1.86E-5; %dynamic viscocity of air (kg/m-s)

%case inputs
T = c.T_C+273.15; %liquid temperature (K)
P = c.P; %operating (headspace) pressure (atm)
Q = c.Q; %liquid flowrate (m3/s)
A_c = c.A_c; %column cross sectional area (m2)
A_orifice = pi*(c.d_orifice_mm/1000/2)^2*c.N_o; %total orifice area (m2)
rho_g = P*MW_Air/(R*T)/1000; %density of gas (kg/m3)
H_CO2 = p.H_CO2*c.H_scale; %Henry's Constants (atm-L/g)
H_O2 = p.H_O2*c.H_scale;
H_N2 = p.H_N2*c.H_scale;

M = length(GQMOM_model_x); %number of quadrature nodes
xspan = 0:0.01:c.towerHeight; %location (m)
v0 = Q/A_orifice; %initial droplet velocity (m/s)
y0 = [repmat([v0,c.c0_CO2,c.c0_O2,c.c0_N2,0],1,M),c.ppm0_CO2,c.ppm0_O2,c.ppm0_N2];

constants = [g, rho_g, p.rho_l, mu_g, p.DAB_CO2, p.DAB_O2, p.DAB_N2, p.s_tension, P, H_CO2, H_O2, H_N2];

%BVP solved with bvp4c
options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NMax', 10000);
if ~isempty(y_guess) && size(y_guess,1) == length(xspan)
    solinit = bvpinit(xspan, @(x) interp1(xspan, y_guess, x)');
else
    solinit = bvpinit(xspan, y0);
end
sol = bvp4c(@(x,y) ODEsys_desorber(x, y, GQMOM_model_x, GQMOM_model_w, constants), ...
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

%liquid holdup
tower_v_tot = A_c*c.towerHeight; %tower volume (m3)
droplet_v_tot_m3 = Q*t_x(end); %total droplet (LH) volume (m3)
liquid_holdup_tot = droplet_v_tot_m3/tower_v_tot; %total liquid holdup (unitless)

%desorption efficiency (liquid side)
desorption_eff_CO2 = (c.c0_CO2-c_CO2_x(end))/c.c0_CO2*100;
desorption_eff_O2 = (c.c0_O2-c_O2_x(end))/c.c0_O2*100;
desorption_eff_N2 = (c.c0_N2-c_N2_x(end))/c.c0_N2*100;

%desorption rates (g/s)
CO2_desorption_rate = (c_CO2_x(1)-c_CO2_x(end))*Q*1000;
O2_desorption_rate = (c_O2_x(1)-c_O2_x(end))*Q*1000;
N2_desorption_rate = (c_N2_x(1)-c_N2_x(end))*Q*1000;

out = struct('caseName', c.caseName, 'stage', 'des', 'sorbent', c.sorbent, 'mdot_liquid', Q*p.rho_l, ...
    'eff_CO2', desorption_eff_CO2, 'eff_O2', desorption_eff_O2, 'eff_N2', desorption_eff_N2, ...
    'CO2_desorption_rate', CO2_desorption_rate, 'O2_desorption_rate', O2_desorption_rate, ...
    'N2_desorption_rate', N2_desorption_rate, 'liquid_holdup_tot', liquid_holdup_tot, ...
    'xspan', xspan(:), 'v_x', v_x, 't_x', t_x, 'c_CO2_x', c_CO2_x, 'c_O2_x', c_O2_x, 'c_N2_x', c_N2_x, ...
    'props', p, 'inputs', c);

end
