function v_eci = ecef2eci(v_ecef, gmst_rad)
    
    C = [cos(gmst_rad)  sin(gmst_rad)  0;
           -sin(gmst_rad) cos(gmst_rad)  0;
           0             0              1]';
    v_eci = C*v_ecef;
end