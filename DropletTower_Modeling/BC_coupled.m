function res = BC_coupled(ya, yb, v0, c0_CO2,c0_O2,c0_N2, ppm0_CO2, ppm0_O2, ppm0_N2, M)

res = zeros(5*M + 3, 1);

for i = 1:M
    idx = (i-1)*5 + 1;
    res(idx) = ya(idx) - v0; % v(0|ri) = v0
    res(idx+1) = ya(idx+1) - c0_CO2; % c_co2(0|ri) = 0
    res(idx+2) = ya(idx+2) - c0_O2; % c_o2(0|ri) = 0
    res(idx+3) = ya(idx+3) - c0_N2; % c_n2(0|ri) = 0
    res(idx+4) = ya(idx+4) - 0; % t(0|ri) = 0
end

res(end-2) = yb(end-2) - ppm0_CO2; % ppm(L) 
res(end-1) = yb(end-1) - ppm0_O2; % ppm(L) 
res(end) = yb(end) - ppm0_N2; % ppm(L) 

end