function p = props_water(T)
%water properties at liquid temperature T (K)

MW_CO2 = 44.01; %CO2 gas molar mass (g/mol)
MW_O2 = 32; %O2 gas molar mass (g/mol)
MW_N2 = 28; %N2 gas molar mass (g/mol)

p.rho_l = 997; %density of water (kg/m3)
p.s_tension = 235.8*(1-T/647.15)^1.256*(1-0.625*(1-T/647.15)); %water surface tension (g/s2)
p.DAB_CO2 = 13.942E-9*(T/227-1)^1.7094; %diffusion coefficient of CO2 in water (m2/s) Lu et al. 2013
p.DAB_O2 = 2.01E-9; %diffusion coefficient of O2 (m2/s)
p.DAB_N2 = 2.00E-9; %diffusion coefficient of N2 (m2/s)

p.H_CO2 = (1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T)))/55342*1000/MW_CO2; %Henry's Constant for CO2 (atm-L/g)
p.H_O2 = (1/(exp(-171.2542+8391.24/T+23.24323*log(T)+0*T)))/55342*1000/MW_O2; %Henry's Constant for O2 (atm-L/g)
p.H_N2 = (1/(exp(-181.587+8632.13/T+24.7981*log(T)+0*T)))/55342*1000/MW_N2; %Henry's Constant for N2 (atm-L/g)

end
