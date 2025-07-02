function gmst_rad = gmst_rad_from_jd(jd)

    % Calculates Greenwich Mean Sidereal Time (GMST) in radians
    % at a given Julian Date (jd).
    % This is a simplified calculation, more precise models exist.

    % J2000 epoch Julian Date
    JD_2000 = 2451545.0;

    % Days since J2000
    d = jd - JD_2000;

    % GMST in seconds (from Astronomical Almanac, simplified)
    gmst_sec = 24110.54841 + 8640184.812866 * d + 0.093104 * (d / 36525)^2 - 6.2e-6 * (d / 36525)^3;

    % Convert to hours and then to radians
    gmst_hours = mod(gmst_sec / 3600, 24); % Modulo 24 hours
    gmst_rad = deg2rad(gmst_hours * 15);   % 15 deg/hour

end