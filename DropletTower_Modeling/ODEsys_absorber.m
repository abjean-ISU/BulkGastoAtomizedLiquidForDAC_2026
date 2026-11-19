function dydx = ODEsys_absorber(x,y,r,w,constants)
%counter-flow absorber, Angelo et al. mass transfer coefficient
%x: position (m)
%y: [v_r1,cCO2_r1,cO2_r1,cN2_r1,t_r1, ... ,v_rM,cCO2_rM,cO2_rM,cN2_rM,t_rM, ppm_CO2, ppm_O2, ppm_N2]
%v: velocity (m/s); c: dissolved concentration (g/L); t: time (s); ppm: gas mole fraction (mole ppm)
%r: droplet radius (m); w: quadrature weights
%constants = [g, rho_g, rho_l, mu_g, DAB_CO2, DAB_O2, DAB_N2, s_tension, P, H_CO2, H_O2, H_N2, A_c, Q, lambda, MW_CO2, MW_O2, MW_N2]

g = constants(1); %acceleration due to gravity (m/s2)
rho_g = constants(2); %density of air (kg/m3)
rho_l = constants(3); %density of liquid (kg/m3)
mu_g = constants(4); %viscocity of air (kg/m-s)
DAB_CO2 = constants(5); %diffusion coefficients (m2/s)
DAB_O2 = constants(6);
DAB_N2 = constants(7);
s_tension = constants(8); %surface tension (g/s2)
P = constants(9); %operating pressure (atm)
H_CO2 = constants(10); %Henry's Constants (atm-L/g)
H_O2 = constants(11);
H_N2 = constants(12);
A_c = constants(13); %column cross sectional area (m2)
Q = constants(14); %liquid flowrate (m3/s)
lambda = constants(15); %1e6/(Ug*P/(R*T)), ppm per (mol/m3/s) per m
MW_CO2 = constants(16); %gas molar masses (g/mol)
MW_O2 = constants(17);
MW_N2 = constants(18);

M = length(r); %number of quadrature points

ppm_CO2 = y(end-2);
ppm_O2 = y(end-1);
ppm_N2 = y(end);

Csat_CO2 = ppm_CO2*P/(10^6*H_CO2); %saturation concentration (g/L)
Csat_O2 = ppm_O2*P/(10^6*H_O2);
Csat_N2 = ppm_N2*P/(10^6*H_N2);

dydx = zeros(5*M+3, 1);

mdot_sum_CO2 = 0;
mdot_sum_O2 = 0;
mdot_sum_N2 = 0;
phi = w(:).*r(:).^3/sum(w(:).*r(:).^3); %liquid volume fraction carried by each size class

for i = 1:M
    idx = (i-1)*5+1;
    v_i = y(idx);
    c_CO2_i = y(idx+1);
    c_O2_i = y(idx+2);
    c_N2_i = y(idx+3);
    r_i = r(i);

    %size class i dynamic quantities
    Re_i = 2*rho_g*v_i*r_i/mu_g; %Reynold's Number
    a_i = 3/r_i; %surface area to volume ratio 1/m
    tau_osc_i = (pi/4)*sqrt(((rho_l*1000)*((2*r_i)^3))/s_tension); %oscillation timescale (s)
    kL_CO2_i = 2*sqrt(DAB_CO2/(pi*tau_osc_i)); %mass transfer coefficient m/s
    kL_O2_i = 2*sqrt(DAB_O2/(pi*tau_osc_i));
    kL_N2_i = 2*sqrt(DAB_N2/(pi*tau_osc_i));

    %liquid phase governing equations
    Gv_i = g - (9/2)*(mu_g/(r_i^2*rho_l))*v_i*(1+0.158*Re_i^(2/3)); %acceleration (m/s2)
    Gc_CO2_i = kL_CO2_i*a_i*(Csat_CO2-c_CO2_i); %rate of change of concentration (g/L/s)
    Gc_O2_i = kL_O2_i*a_i*(Csat_O2-c_O2_i);
    Gc_N2_i = kL_N2_i*a_i*(Csat_N2-c_N2_i);

    dydx(idx) = Gv_i/v_i;
    dydx(idx+1) = Gc_CO2_i/v_i;
    dydx(idx+2) = Gc_O2_i/v_i;
    dydx(idx+3) = Gc_N2_i/v_i;
    dydx(idx+4) = 1/v_i;

    %volume-weighted gas phase source, each class at its own velocity
    mdot_sum_CO2 = mdot_sum_CO2 + phi(i)*Gc_CO2_i/v_i;
    mdot_sum_O2  = mdot_sum_O2  + phi(i)*Gc_O2_i/v_i;
    mdot_sum_N2  = mdot_sum_N2  + phi(i)*Gc_N2_i/v_i;
end

%gas phase governing equations
mdot_CO2 = Q/A_c*mdot_sum_CO2; %g/L/s
mdot_O2  = Q/A_c*mdot_sum_O2;
mdot_N2  = Q/A_c*mdot_sum_N2;
dydx(end-2) = lambda*mdot_CO2*1000/MW_CO2; %g/L/s -> mol/m3/s -> ppm/m
dydx(end-1) = lambda*mdot_O2*1000/MW_O2;
dydx(end)   = lambda*mdot_N2*1000/MW_N2;

end
