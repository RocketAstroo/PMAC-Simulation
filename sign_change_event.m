function [value, isterminal, direction] = sign_change_event(t, state, p)

    r = state(1:3);      
    v = state(4:6);    
    e = state(7:10); 
    e = e/norm(e); %Normalize Quaternion to avoid drifiting due to accumulation of errors
    w = state(11:13);   

    %Translational dynamics
    rcap = r/norm(r);
    force_gravity = ((p.gravitational_const*p.mass_earth*p.mass_sat)/(norm(r)^2))*(-rcap);
    accel = force_gravity/p.mass_sat;

    %Rotational dynamics
    edot = (1/2)*[0, w(3), -w(2), w(1); ...
                 -w(3), 0, w(1), w(2); ...
                  w(2), -w(1), 0, w(3); ...
                  -w(1), -w(2), -w(3), 0]*e;
    
    %Calling IGRF model
    current_jd = p.jd + t/(24*3600);
    gmst_rad = gmst_rad_from_jd(current_jd);
    
    r_ecef = eci2ecef(r, gmst_rad); 
    lat_long = ecef2lat_long(r_ecef);

    [Bn_nT, Be_nT, Bd_nT] = igrf("1-Jul-2025", lat_long(1)*(180/pi), lat_long(2)*(180/pi), (p.semi_major_orbit)/1000,'geocentric');
    Bu_nT = -Bd_nT;
    Bn = Bn_nT*1e-9;
    Be = Be_nT*1e-9;
    Bu = Bu_nT*1e-9;
    B_neu = [Bn; Be; Bu];
    theta3 = 90 + atan(r_ecef(2)/r_ecef(1));
    theta1 = 90 - atan(r_ecef(3)/(sqrt(r_ecef(1)^2 + r_ecef(2)^2)));
    B_ecef = neu2ecef(B_neu,theta3,theta1);
    B_eci = ecef2eci(B_ecef, gmst_rad);

    psi = atan(2*(e(1)*e(2) + e(4)*e(3))/(1-2*(e(2)^2 + e(3)^2)));
    theta = asin(2*(e(4)*e(2) - e(1)*e(3)));
    phi = atan(2*(e(2)*e(3) + e(4)*e(1))/(1-2*(e(1)^2 + e(2)^2)));

    B_body = eci2body(B_eci, psi, theta, phi);
    H_body = B_body/p.mu0;

    fprintf('Time: %f s, H_body magnitude: %e A/m\n', t, norm(H_body));

    dH_dt_vector = -cross(w, H_body); 

     dH_dt_epsilon = 1e-9;
     value = zeros(3,1);

      B_hyst = (2/pi)*p.Bs*atan(p.p0*(H_body + p.Hc* p.current_Hc_signs));

    m_hyst = B_hyst*(p.rod_vol/p.mu0);
    m_hyst_tot = m_hyst.*p.no_rods;
    moment_hyst = cross(m_hyst_tot, B_body);

    if p.current_Hc_signs(1) == 1 % Current Hc is +Hc, means dH/dt was < 0. Looking for dH/dt to cross +epsilon (i.e., dH/dt becomes > 0)
        value(1) = dH_dt_vector(1) - dH_dt_epsilon; % Event when dH_dt_vector(1) = dH_dt_epsilon
    else % p.current_Hc_signs(1) == -1. Current Hc is -Hc, means dH/dt was > 0. Looking for dH/dt to cross -epsilon (i.e., dH/dt becomes < 0)
        value(1) = dH_dt_vector(1) + dH_dt_epsilon; % Event when dH_dt_vector(1) = -dH_dt_epsilon
    end

    % Event for Y-component
    if p.current_Hc_signs(2) == 1 
        value(2) = dH_dt_vector(2) - dH_dt_epsilon; 
    else 
        value(2) = dH_dt_vector(2) + dH_dt_epsilon; 
    end

    % Event for Z-component
    if p.current_Hc_signs(3) == 1 
        value(3) = dH_dt_vector(3) - dH_dt_epsilon; 
    else 
        value(3) = dH_dt_vector(3) + dH_dt_epsilon; 
    end


    
    % When an event occurs, terminate the integration segment
    isterminal = [1; 1; 1]; % 1 = terminate integration for the detected events (Y and Z rods)
                         % This array must match the size of 'value'
    % Detect zero-crossings from any direction (H can increase then decrease, or vice-versa)
    direction = [0; 0; 0];  % 0 = detect all zero-crossings (up or down)
                         % This array must match the size of 'value'

      % Optional: Debugging fprintf statements for values leading to event
    fprintf('Event Func: t=%.10f, dH_x/dt=%.4e, dH_y/dt=%.4e, dH_z/dt=%.4e, Hc_signs=[%.0f,%.0f, %.0f], EventValue=[%.4e, %.4e, %.4e]\n', ...
     t, dH_dt_vector(1), dH_dt_vector(2), dH_dt_vector(3), p.current_Hc_signs(1), p.current_Hc_signs(2), p.current_Hc_signs(3), ...
     value(1), value(2), value(3));

   

     power_hyst_dissipated = dot(moment_hyst, w);
     fprintf('Time: %f s, Power_Hyst: %e W\n', t, power_hyst_dissipated);




end