function jd = gregorian_to_julian(date_vec)
    % Simplified Julian Date calculation
    % [year, month, day, hour, minute, second]
    year = date_vec(1); month = date_vec(2); day = date_vec(3);
    hour = date_vec(4); minute = date_vec(5); second = date_vec(6);

    if month <= 2
        year = year - 1;
        month = month + 12;
    end

    A = floor(year / 100);
    B = 2 - A + floor(A / 4);

    jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;
    jd = jd + (hour + minute/60 + second/3600) / 24;
end