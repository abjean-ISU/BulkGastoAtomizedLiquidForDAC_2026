function [GQMOM_model_w, GQMOM_model_x] = dropletQuadrature(DropletSize)
%weights and abscissae (radii, m) for the droplet size distribution
%DropletSize: measured droplet diameters (mm) OR a single mean diameter (mm)

if isscalar(DropletSize)
    GQMOM_model_x = DropletSize/1000/2; %convert mean drop diameter to radius (m)
    GQMOM_model_w = 1; %using only 1 abscissae
else
    %normalizing diameters against maximum diameter
    diameter_max = max(DropletSize);
    diameter_norm = DropletSize./diameter_max;

    %computing the histogram and moments
    [hist_weights,edges] = histcounts(diameter_norm,20,'Normalization','probability');
    hist_abscissae = edges(1:end-1)+ diff(edges)/2;
    for i = 0:10
        calc = hist_weights.*(hist_abscissae.^i);
        hist_moments(i+1) = sum(calc);
    end

    %4-node quadrature from moments 0-7
    [GQMOM_model_w,GQMOM_model_x,GQMOM_model_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),4);
    if any(GQMOM_model_werror ~= 0)
        warning('dropletQuadrature: PD_gqmom_0_1 returned werror ~= 0')
    end

    %converting x to radii in m
    GQMOM_model_x = GQMOM_model_x*diameter_max/2/1000; %m
end

end
