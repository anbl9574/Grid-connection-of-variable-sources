function tiltvar(SITE, TILTS, AZIMUTH, ALBEDO)

% Monthly time-keeping
days = [1 31 28 31 30 31 30 31 31 30 31 30 31];
days_cum = cumsum(days);

% Summer months
ind1 = (days_cum(4)-1)*24 + 1;
ind2 = (days_cum(10)-1)*24;

% Summer and winter indices
inds = ind1:ind2;
indw = setdiff(1:8760, inds);


ind = 1;

for i = TILTS
    
    [G1 G2 G3] = solrad(SITE, i, AZIMUTH, ALBEDO, false);
    
    G = (G1 + G2 + G3) / 1000;
    
    Gs(ind) = sum(G(inds));
    Gw(ind) = sum(G(indw));  
    Gy(ind) = sum(G(1:8760));
    
    ind = ind + 1;
    
end

figure
p = plot(TILTS, Gs, '-ro', TILTS, Gw, '-bo', TILTS, Gy, '-ko');
xlabel('Tilt angle (\circ)')
ylabel('kWh/m^2')
legend({'Apr-Sep', 'Oct-Mar', 'Whole year'})
title(['Total in-plane irradiance at azimuth ' num2str(round(AZIMUTH)) '\circ as function of tilt angle'])
set(p(1), 'LineWidth', 2)
set(p(2), 'LineWidth', 2)
set(p(3), 'LineWidth', 2)

