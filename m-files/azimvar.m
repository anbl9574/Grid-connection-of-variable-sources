function azimvar(SITE, TILT, AZIMUTHS, ALBEDO)

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

for i = AZIMUTHS
    
    [G1 G2 G3] = solrad(SITE, TILT, i, ALBEDO, false);
    
    G = (G1 + G2 + G3) / 1000;
    
    Gs(ind) = sum(G(inds));
    Gw(ind) = sum(G(indw));  
    Gy(ind) = sum(G(1:8760));
    
    ind = ind + 1;
    
end

figure
p = plot(AZIMUTHS, Gs, '-ro', AZIMUTHS, Gw, '-bo', AZIMUTHS, Gy, '-ko');
xlabel('Azimuth angle (\circ)')
ylabel('kWh/m^2')
legend({'Apr-Sep', 'Oct-Mar', 'Whole year'})
title(['Total in-plane irradiance at tilt ' num2str(round(TILT)) '\circ as function of azimuth angle'])
set(p(1), 'LineWidth', 2)
set(p(2), 'LineWidth', 2)
set(p(3), 'LineWidth', 2)