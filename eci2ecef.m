function v_ecef = eci2ecef(v_eci, gmst_rad)
    
    C = [cos(gmst_rad)  sin(gmst_rad)  0;
           -sin(gmst_rad) cos(gmst_rad)  0;
           0             0              1];
    v_ecef = C*v_eci;
end