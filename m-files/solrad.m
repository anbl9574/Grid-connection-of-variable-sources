function [beamTilted diffuseTilted groundrefTilted] = solrad(SITE, PANEL_TILT, PANEL_AZIMUTH, ALBEDO, PLOT)

% Location constants
LATITUDE  = SITE.latitude;
LONGITUDE = SITE.longitude;
TIME_ZONE = SITE.timeZone;

% Radiation components
globalHor  = SITE.global;
diffuseHor = SITE.diffuse;
beamHor    = globalHor - diffuseHor;

% Initialize vectors for new radiation components
extraterr       = zeros(8760,1);
beamTilted      = zeros(8760,1);
diffuseTilted   = zeros(8760,1);
groundrefTilted = zeros(8760,1);

% Time
hour      = SITE.hourOfDay;
dayOfYear = SITE.dayOfYear;

% Time step in minutes
TIME_STEP = 60; 

% Calculate radiation on tilted planes
% for all time steps of the year
for i = 1:8760 

    % Minute of day in standard time (midpoint of hour)
    minOfDay = TIME_STEP * (hour(i) - 0.5);
    
    % Standard meridian
    longitudeStd = -TIME_ZONE * 15;
    
    
    % ### SOLAR TIME ###
    
    % Equation of time
    B = 2 * pi * (dayOfYear(i) - 1) / 365;
    E = 229.18 * ( 0.000075 + 0.001868*cos(B) - 0.032077*sin(B) - 0.014615*cos(2*B) - 0.04089*sin(2*B) ); 

    % Correction for longitude and equation of time
    minOfDaySolar = minOfDay - 4 * (longitudeStd - LONGITUDE) + E;
    
    
    % ### BEAM RADIATION ###
    
    % (Re)define all parameters in radians
    lat   = LATITUDE * 2 * pi / 360;                                    % Latitude
    decl  = 2 * pi * 23.45 * sin(2*pi*(284 + dayOfYear(i))/365) / 360;  % Declination
    tilt  = 2 * pi * PANEL_TILT / 360;                                  % Panel tilt
    azim  = 2 * pi * PANEL_AZIMUTH / 360;                               % Panel azimuth
    hang  = pi/12 * (minOfDaySolar/60 - 12);                            % Hour angle
    
    % Calculate cosine of angle of incidence on tilted plane
    cosTheta = sin(decl).*sin(lat).*cos(tilt) - sin(decl).*cos(lat).*sin(tilt).*cos(azim) + cos(decl).*cos(lat).*cos(tilt).*cos(hang) + cos(decl).*sin(lat).*sin(tilt).*cos(azim).*cos(hang) + cos(decl).*sin(tilt).*sin(azim).*sin(hang);

    % Calculate cosine of angle of incidence on horizontal plane
    cosThetaZenith = cos(decl).*cos(lat).*cos(hang) + sin(decl).*sin(lat);
    
    
    % ### GEOMETRIC FACTOR ###
    
    if cosTheta > 0 && cosThetaZenith > 0
        Rb = cosTheta / cosThetaZenith;
    else
        Rb = 0;
    end
    
    
    % ### EXTRATERRESTRIAL RADIATION AND ANISOTROPY INDEX ###
    
    extraterr(i) = 1367 * (1 + 0.033*cos(2*pi*dayOfYear(i)/365)) * (cos(lat)*cos(decl)*cos(hang) + sin(lat)*sin(decl));
  
    if extraterr(i) > 0
        Ai = beamHor(i) / extraterr(i);
    else
        Ai = 0;
    end
    
    
    % ###  RADIATION COMPONENTS ON THE TILTED PLANE  ###

    % Beam radiation
    beamTilted(i) = Rb * beamHor(i);

    % Diffuse radiation
    diffuseTilted(i) = ((1 - Ai) * (1 + cos(tilt))/2 + Ai * Rb) * diffuseHor(i);

    % Ground-reflected radiation
    groundrefTilted(i) = ALBEDO * (1 - cos(tilt))/2 * (beamHor(i) + diffuseHor(i));
 
    % Global tilted radiation
    globalTilted = beamTilted + diffuseTilted + groundrefTilted;
    
end


if PLOT
    
    % Monthly time-keeping
    days = [1 31 28 31 30 31 30 31 31 30 31 30 31];
    days_cum = cumsum(days);
    
    % Summer months
    ind1 = (days_cum(4)-1)*24 + 1;
    ind2 = (days_cum(10)-1)*24;
 
    % Summer and winter indices
    inds = ind1:ind2;
    indw = setdiff(1:8760, inds);
    
    % Daily mean radiation curves
    mrs = mean_load_curve(globalTilted(inds)', 60);
    mrw = mean_load_curve(globalTilted(indw)', 60);
    
    % Monthly totals
    for i = 1:12
        ind1 = (days_cum(i)-1)*24 + 1;
        ind2 = (days_cum(i+1)-1)*24;
        sm(i,1) = sum(globalTilted(ind1:ind2))/1000;
    end

    f = figure;
    
    % Plot hourly values
    h1 = subplot(3,1,1); plot(1:8760, globalTilted, '-k')
    xlim([1 8760])
    ylim([0 1200])
    set(h1, 'XTick', [1, 8760])
    xlabel('Hour of year')
    ylabel('W/m^2')
    title('Hourly global radiation on the tilted plane')
    annotation(f,'textbox',[0 0 1 1],'FitBoxToText','on','EdgeColor','none','String',['Total in-plane irradiance: ' num2str(round(sum(globalTilted/1000))) ' kWh/m2/year']);
    
    % Plot mean daily patterns
    h2 = subplot(3,1,2); p = plot(1:24, mrs, '-ko', 1:24, mrw, '-kx');
    xlim([1 24])
    ylim([0 1000])
    set(h2, 'XTick', [1, 6, 12, 18, 24])
    xlabel('Hour of day')
    ylabel('W/m^2')
    legend({'Apr-Sep', 'Oct-Mar'})
    title('Hourly averages of global radiation on the tilted plane')
    set(p(1), 'LineWidth', 2)
    set(p(2), 'LineWidth', 2)
    
    % Plot monthly totals
    h3 = subplot(3,1,3); p = plot(1:12, sm, '-kx');
    xlim([0 13])
    ylim([0 300])
    set(h3, 'XTick', 1:12)
    xlabel('Month')
    ylabel('kWh/m^2/month')
    set(h3, 'XTickLabel', {'J','F','M','A','M','J','J','A','S','O','N','D'})
    title('Monthly totals of global radiation on the tilted plane')
    set(p, 'LineWidth', 2)
    
end