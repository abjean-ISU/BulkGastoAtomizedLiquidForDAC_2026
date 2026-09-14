function res = BC_coupled_desorption(ya, yb, v0, c0_CO2, c0_O2, c0_N2, ppm_vac, M)

res = zeros(5*M + 3, 1);

for i = 1:M
    idx = (i-1)*5 + 1;
    res(idx) = ya(idx) - v0; % v(0|ri) = v0
    res(idx+1) = ya(idx+1) - c0_CO2; % c_co2(0|ri) = 0
    res(idx+2) = ya(idx+2) - c0_O2; % c_o2(0|ri) = 0
    res(idx+3) = ya(idx+3) - c0_N2; % c_n2(0|ri) = 0
    res(idx+4) = ya(idx+4) - 0; % t(0|ri) = 0
end

% Fix ppm at BOTH ends to ppm_vac
% This represents vacuum pump maintaining near-zero ppm throughout
% Gas inlet (bottom, x=L):
res(end-2) = yb(end-2) - ppm_vac;  % ppm_CO2(L) = ppm_vac
res(end-1) = yb(end-1) - ppm_vac;  % ppm_O2(L)  = ppm_vac
res(end)   = yb(end)   - ppm_vac;  % ppm_N2(L)  = ppm_vac

end