function p = props_IL(IL, T)
%ionic liquid properties at liquid temperature T (K)
%IL: 'bmimTf2N', 'emimTf2N', 'bmimDCA', 'bmimSCN'

MW_CO2 = 44.01; %CO2 gas molar mass (g/mol)
MW_O2 = 32; %O2 gas molar mass (g/mol)
MW_N2 = 28; %N2 gas molar mass (g/mol)
R_J = 8.314; %ideal gas constant (J/mol/K)
T_ref = 298.15; %K Reference temperatures for Arrhenius equations

switch IL
    case 'bmimTf2N'
        MW_IL = 419.36; %g/mol [bmim][Tf2N]
        rho_l = -1.89*T+1978; %IL density (kg/m3) DOI: 10.1016/j.ijggc.2016.08.022
        s_tension = -0.039*T+44.77; %surface tension (g/s2) DOI: 10.1016/j.jct.2018.12.036
        eta_l = 45*exp(32000/R_J*(1/T-1/T_ref)); %dynamic viscocity (cP) DOI: 10.1021/ie070501u

        %Henry's Constant for CO2 DOI: 10.1021/ie070501u
        rho_l_CO2ref = -1.89*T_ref+1978; %IL density (kg/m3) DOI: 10.1016/j.ijggc.2016.08.022
        Href_CO2 = 34*0.98*MW_IL/rho_l_CO2ref/MW_CO2; %reference H converted from bar to atm-L/g
        H_CO2 = Href_CO2*exp(-12000/R_J*(1/T-1/T_ref)); %Henry's Constant (atm-L/g)

        %Henry's Constant for O2 DOI: 10.1021/jp046404l
        rho_l_O2ref = -1.89*323+1978; %IL density (kg/m3) DOI: 10.1016/j.ijggc.2016.08.022
        Href_O2 = 1110*0.98*MW_IL/rho_l_O2ref/MW_O2; %reference H converted from MPa to atm-L/g
        H_O2 = Href_O2*exp(20600/R_J*(1/T-1/323)); %Henry's Constant (atm-L/g)

        %Henry's Constant for N2 DOI: 10.1021/acs.iecr.5c03061
        rho_l_N2ref = -1.89*308+1978; %IL density (kg/m3) DOI: 10.1016/j.ijggc.2016.08.022
        Href_N2 = 700*0.98*MW_IL/rho_l_N2ref/MW_N2; %reference H converted from bar to atm-L/g
        H_N2 = Href_N2; %enthalpy of absorption assumed 0 according to 10.1021/acs.iecr.5c03061 and 10.1021/acs.jced.3c00159

    case 'emimTf2N'
        MW_IL = 391.31; %g/mol [emim][Tf2N]
        rho_l = -0.89253*T+1814.3; %IL density (kg/m3) DOI: 10.1016/j.fluid.2021.113211
        s_tension = -0.12*T+77.49; %surface tension (g/s2) DOI: 10.1016/j.jct.2013.02.021
        eta_l = 0.12*exp(892.53/(T-139.12)); %dynamic viscocity (cP) DOI: 10.1016/j.fluid.2021.113211

        %Henry's Constant for CO2 DOI: 10.1021/ie0704142
        rho_l_CO2ref = -0.9452*298.15+1802; %IL density (kg/m3) DOI: 10.1021/je034261a
        Href_CO2 = 39*MW_IL/rho_l_CO2ref/MW_CO2; %reference H converted from atm to atm-L/g
        H_CO2 = Href_CO2*exp(-13000/R_J*(1/T-1/298.15)); %Henry's Constant (atm-L/g)

        %Henry's Constant for O2 DOI: 10.1021/jp046404l - same as [bmim][Tf2N]
        rho_l_O2ref = -0.9452*323+1802; %IL density (kg/m3) DOI: 10.1021/je034261a
        Href_O2 = 1110*0.98*MW_IL/rho_l_O2ref/MW_O2; %reference H converted from MPa to atm-L/g
        H_O2 = Href_O2*exp(20600/R_J*(1/T-1/323)); %Henry's Constant (atm-L/g)

        %Henry's Constant for N2 DOI: 10.1021/acs.jced.3c00159
        rho_l_N2ref = -0.9452*333.2+1802; %IL density (kg/m3) DOI: 10.1021/je034261a
        Href_N2 = 1090*MW_IL/rho_l_N2ref/MW_N2; %reference H converted from atm to atm-L/g
        H_N2 = Href_N2*exp(500/R_J*(1/T-1/333.2)); %Henry's Constant (atm-L/g)

    case 'bmimDCA'
        MW_IL = 205.26; %g/mol [bmim][DCA]
        rho_l = -0.55*T+1223; %IL density (kg/m3) DOI: 10.1021/je034261a
        s_tension = (-0.0775E-3*T+71.88E-3)*1000; %surface tension (g/s2) DOI: 10.6100/IR639177
        eta_l = 0.05417*exp(1000/(T-138)); %dynamic viscocity (cP) DOI: 10.6100/IR639177

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
        H_N2 = Href_N2; %enthalpy of absorption assumed 0 according to 10.1021/acs.iecr.5c03061 and 10.1021/acs.jced.3c00159

    case 'bmimSCN'
        MW_IL = 197.3; %g/mol [bmim][SCN]
        rho_l = (-5.8E-4*T+1.2426)*1000; %IL density (kg/m3) DOI: 10.6100/IR639177
        s_tension = (-0.2023E-3*T+106.11E-3)*1000; %surface tension (g/s2) DOI: 10.6100/IR639177
        eta_l = -.8075*T+285.97; %dynamic viscocity (cP) DOI: 10.1021/je8008254

        %Henry's Constant for CO2 DOI: 10.6100/IR639177
        rho_l_CO2ref = (-5.8E-4*303+1.2426)*1000; %IL density (kg/m3) DOI: 10.6100/IR639177
        Href_CO2 = 101*0.98*MW_IL/rho_l_CO2ref/MW_CO2; %reference H converted from bar to atm-L/g
        H_CO2 = Href_CO2*exp(-10170/R_J*(1/T-1/303)); %Henry's Constant (atm-L/g)

        %Henry's Constant for O2 DOI: 10.1021/jp046404l - same as [bmim][Tf2N]
        rho_l_O2ref = (-5.8E-4*323+1.2426)*1000; %IL density (kg/m3) DOI: 10.6100/IR639177
        Href_O2 = 1110*0.98*MW_IL/rho_l_O2ref/MW_O2; %reference H converted from MPa to atm-L/g
        H_O2 = Href_O2*exp(20600/R_J*(1/T-1/323)); %Henry's Constant (atm-L/g)

        %Henry's Constant for N2 DOI: 10.1021/acs.iecr.5c03061
        rho_l_N2ref = (-5.8E-4*308+1.2426)*1000; %IL density (kg/m3) DOI: 10.6100/IR639177
        Href_N2 = 9800*0.98*MW_IL/rho_l_N2ref/MW_N2; %reference H converted from bar to atm-L/g
        H_N2 = Href_N2; %enthalpy of absorption assumed 0 according to 10.1021/acs.iecr.5c03061 and 10.1021/acs.jced.3c00159

    otherwise
        error('props_IL: unknown ionic liquid %s', IL)
end

DAB_CO2 = (6.7E5*(eta_l^-0.66)*(MW_IL^-0.89)*(rho_l/1000)^4.8*T^-3.3)/10000; %CO2 diffusivity (m2/s) DOI: 10.1021/ie070501u

p.rho_l = rho_l;
p.s_tension = s_tension;
p.DAB_CO2 = DAB_CO2;
p.DAB_O2 = DAB_CO2; %approximation
p.DAB_N2 = DAB_CO2; %approximation
p.H_CO2 = H_CO2;
p.H_O2 = H_O2;
p.H_N2 = H_N2;

end
