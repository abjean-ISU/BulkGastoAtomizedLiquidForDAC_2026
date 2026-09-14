function dydx = ODEsys_Angelo_Desorption(x,y,r,w,constants)
%inputs
%x: position 
%y[v_r1,c_r1,t_r1,v_r2,c_r2,t_r2,...v_rM,c_rM,t_rM,ppm]
%v: velocity (m/s) v(0)=Q/A
%c: concentration (g/L) c(0)=0
%t: time (s) t(0)=0
%ppm: concentration of CO2 in gas strean (ppm)
%r: droplet radius (m)
%w: quadrature weights 
%constants = [g, rho_g, rho_l, mu_g, DAB_CO2, DAB_O2, DAB_N2, s_tension, P, H_CO2, H_O2, H_N2, A_c, Q, lambda]; 

%constants 
g = constants(1); %acceleration due to gravity (m/s2)
rho_g = constants(2); %density of air (kg/m3)
rho_l = constants(3); %density of water (kg/m3) 
mu_g = constants(4); %viscocity of air (m2/s)
DAB_CO2 = constants(5); %diffusion coefficient of CO2 in water at 25C (m2/s)
DAB_O2 = constants(6); 
DAB_N2 = constants(7); 
s_tension = constants(8); %water surface tenstion (g/s2)
P = constants(9); %operating pressue (atm)
H_CO2 = constants(10); %Henry's Constant (atm-L/g)
H_O2 = constants(11);
H_N2 = constants(12);
A_c = constants(13); %Cross Sectional column area (m2)
Q = constants(14); %liquid flowrate (m3/s)
ppm_vac = constants(15); %fixed headspace ppm

M = length(r); %number of quadrature points

ppm_CO2 = y(end-2); %ppm vector
ppm_O2 = y(end-1); %ppm vector
ppm_N2 = y(end); %ppm vector

Csat_CO2 = ppm_CO2*P/(10^6*H_CO2); %saturation concentration (g/L)
Csat_O2 = ppm_O2*P/(10^6*H_O2); %saturation concentration (g/L)
Csat_N2 = ppm_N2*P/(10^6*H_N2); %saturation concentration (g/L)

%initialized ouput vectors
dydx = zeros(5*M+3, 1); 

%initialized mass transfer rate accumulator
mdot_sum_CO2 = 0;
mdot_sum_O2 = 0;
mdot_sum_N2 = 0;

for i = 1:M
    %indexing state vectors
    idx = (i-1)*5+1;
    v_i = y(idx);
    c_CO2_i = y(idx+1);
    c_O2_i = y(idx+2);
    c_N2_i = y(idx+3);
    t_i = y(idx+4);
    r_i = r(i);
    w_i = w(i); 


    %size class i dynamic quantities 
    Re_i = 2*rho_g*v_i*r_i/mu_g; %Reynold's Number 
    a_i = 3/r_i; %surface area to volume ratio 1/m
    tau_osc_i = (pi/4)*sqrt(((rho_l*1000)*((2*r_i)^3))/s_tension); %oscillation timescale (s)
    kL_CO2_i = 2*sqrt(DAB_CO2/(pi*tau_osc_i)); %mass transfer coefficient m/s
    kL_O2_i = 2*sqrt(DAB_O2/(pi*tau_osc_i)); %mass transfer coefficient m/s
    kL_N2_i = 2*sqrt(DAB_N2/(pi*tau_osc_i)); %mass transfer coefficient m/s

    %liquid phase size class governing equations 
    Gv_i = g - (9/2)*(mu_g/(r_i^2*rho_l))*v_i*(1+0.158*Re_i^(2/3)); %acceleration (m/s2)
    Gc_CO2_i = kL_CO2_i*a_i*(Csat_CO2-c_CO2_i); %rate of change of CO2 concentration (g/L/s)
    Gc_O2_i = kL_O2_i*a_i*(Csat_O2-c_O2_i); %rate of change of O2 concentration (g/L/s)
    Gc_N2_i = kL_N2_i*a_i*(Csat_N2-c_N2_i); %rate of change of N2 concentration (g/L/s)

    %liquid phase ODEs for size class
    dydx(idx) = Gv_i/v_i; % Update velocity derivative
    dydx(idx+1) = Gc_CO2_i/v_i; % Update concentration derivative
    dydx(idx+2) = Gc_O2_i/v_i; % Update concentration derivative
    dydx(idx+3) = Gc_N2_i/v_i; % Update concentration derivative
    dydx(idx+4) = 1/v_i; % Update time derivative

    %gas phase governing equation for size class
    mdot_sum_CO2 = mdot_sum_CO2 + w_i*kL_CO2_i*a_i*(Csat_CO2 - c_CO2_i);
    mdot_sum_O2 = mdot_sum_O2 + w_i*kL_O2_i*a_i*(Csat_O2 - c_O2_i); 
    mdot_sum_N2 = mdot_sum_N2 + w_i*kL_N2_i*a_i*(Csat_N2 - c_N2_i); 
end

%local averaged quantites 
droplet_v = sum(w.*(4/3)*pi.*r.^3); %average droplet volume (m3)
v_bar = sum(w.*y(1:5:5*M)); %average droplet velcotiy (m/s)
flux_const = Q/droplet_v; %droplets/s
N_x = flux_const/v_bar; %droplet count droplets/m
n_x = N_x/A_c; %droplet density droplets/m3
epsilon_L = n_x*droplet_v;

%averaged gas phase governing equation
mdot_CO2 = epsilon_L*mdot_sum_CO2; %g/L/s
mdot_O2 = epsilon_L*mdot_sum_O2; %g/L/s
mdot_N2 = epsilon_L*mdot_sum_N2; %g/L/s
dydx(end-2) = 0; %gas-phase concentration maintained at zero (ppm/m)
dydx(end-1) = 0; %gas-phase concentration (ppm/m)
dydx(end) = 0; %gas-phase concentration (ppm/m)

end


