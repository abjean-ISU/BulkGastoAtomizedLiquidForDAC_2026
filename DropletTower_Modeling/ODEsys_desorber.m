function dydx = ODEsys_desorber(x,y,r,w,constants)
%vacuum desorber, Angelo et al. mass transfer coefficient
%headspace gas held at ppm_vac (gas-phase ODEs = 0), liquid-controlled desorption
%y: [v_r1,cCO2_r1,cO2_r1,cN2_r1,t_r1, ... ,t_rM, ppm_CO2, ppm_O2, ppm_N2]
%constants = [g, rho_g, rho_l, mu_g, DAB_CO2, DAB_O2, DAB_N2, s_tension, P, H_CO2, H_O2, H_N2]

g = constants(1); %acceleration due to gravity (m/s2)
rho_g = constants(2); %density of gas (kg/m3)
rho_l = constants(3); %density of liquid (kg/m3)
mu_g = constants(4); %viscocity of gas (kg/m-s)
DAB_CO2 = constants(5); %diffusion coefficients (m2/s)
DAB_O2 = constants(6);
DAB_N2 = constants(7);
s_tension = constants(8); %surface tension (g/s2)
P = constants(9); %operating pressure (atm)
H_CO2 = constants(10); %Henry's Constants (atm-L/g)
H_O2 = constants(11);
H_N2 = constants(12);

M = length(r); %number of quadrature points

Csat_CO2 = y(end-2)*P/(10^6*H_CO2); %saturation concentration (g/L)
Csat_O2 = y(end-1)*P/(10^6*H_O2);
Csat_N2 = y(end)*P/(10^6*H_N2);

dydx = zeros(5*M+3, 1);

for i = 1:M
    idx = (i-1)*5+1;
    v_i = y(idx);
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
    dydx(idx) = Gv_i/v_i;
    dydx(idx+1) = kL_CO2_i*a_i*(Csat_CO2-y(idx+1))/v_i;
    dydx(idx+2) = kL_O2_i*a_i*(Csat_O2-y(idx+2))/v_i;
    dydx(idx+3) = kL_N2_i*a_i*(Csat_N2-y(idx+3))/v_i;
    dydx(idx+4) = 1/v_i;
end

%gas phase held at ppm_vac
dydx(end-2:end) = 0;

end
