%Main file where the simulation actual run

% Clear workspace, close figures, clear command window
clear; close all; clc;

%Getting all the parameters
p = parameters();

%Initialization of the state vector
 %Initialization of the translational parameters
  x0_eci = p.semi_major_orbit;
  y0_eci = 0;
  z0_eci = 0;
  r0_eci = [x0_eci;y0_eci;z0_eci];

  xdot0_eci = 0;
  ydot0_eci = p.vel_orbit*cos(p.inclination_orbit);
  zdot0_eci = p.vel_orbit*sin(p.inclination_orbit);
  v0_eci = [xdot0_eci;ydot0_eci;zdot0_eci];

 %Initialization of the rotational parameters
  %Initialization of euler angles 3-2-1 sequence (intrinsic rotations)
  psi0 = 0; 
  theta0 = 0; 
  phi0 = 0; 

  %Conversion of euler angles to euler parameters/quaternions
  e10 = sin(phi0/2)*cos(theta0/2)*cos(psi0/2) - cos(phi0/2)*sin(theta0/2)*sin(psi0/2);
  e20 = cos(phi0/2)*sin(theta0/2)*cos(psi0/2) + sin(phi0/2)*cos(theta0/2)*sin(psi0/2);
  e30 = cos(phi0/2)*cos(theta0/2)*sin(psi0/2) - sin(phi0/2)*sin(theta0/2)*cos(psi0/2);
  eta0 = cos(phi0/2)*cos(theta0/2)*cos(psi0/2) + sin(phi0/2)*sin(theta0/2)*sin(psi0/2);
  e0 = [e10; e20; e30; eta0];
  e0 = e0/ norm(e0); %For the sake of Normalization of the parameters

  %Initialization of angular velocity of the body w.r.t eci but represented in body axes
  wx0 = 0.07;
  wy0 = 0.07;
  wz0 = 0.07;
  w0 = [wx0; wy0; wz0];
  
  %Initialization of the overall state vector
  state0 = [r0_eci;v0_eci;e0;w0];

%Initialization of date, time and timespan of the simulation
start_jd = gregorian_to_julian(p.date); % Julian Date of simulation start
no_of_orbits = 175;
tspan = [0, 100000];

% Initialize solution storage for all segments
all_t_sol = [];
all_state_sol = [];

% Current state for the segment
t_current_segment_start = tspan(1);
current_state_restart = state0;
signs = [];
last = 0;
prev_last = 0;

%Initialize Hc signs
current_jd = p.jd + tspan(1)/(24*3600);
gmst_rad = gmst_rad_from_jd(current_jd);
    
r_ecef = eci2ecef(state0(1:3), gmst_rad); 
lat_long = ecef2lat_long(r_ecef);
[Bn_nT, Be_nT, Bd_nT] = igrf("1-Jul-2025", lat_long(1)*(180/pi), lat_long(2)*(180/pi), (p.semi_major_orbit)/1000,'geocentric');
        Bu_nT = -Bd_nT;
        Bn = Bn_nT*1e-9;
        Be = Be_nT*1e-9;
        Bu = Bu_nT*1e-9;
        B_neu_ini = [Bn; Be; Bu];
        theta3 = 90 + atan2(r_ecef(2),r_ecef(1));
        theta1 = 90 - atan2(r_ecef(3),(sqrt(r_ecef(1)^2 + r_ecef(2)^2)));
        B_ecef_ini = neu2ecef(B_neu_ini,theta3,theta1);
        B_eci_ini = ecef2eci(B_ecef_ini, gmst_rad);

        e_ini = state0(7:10);
        psi_ini = atan2(2*(e_ini(1)*e_ini(2) + e_ini(4)*e_ini(3)),(1-2*(e_ini(2)^2 + e_ini(3)^2)));
        theta_ini = asin(2*(e_ini(4)*e_ini(2) - e_ini(1)*e_ini(3)));
        phi_ini = atan2(2*(e_ini(2)*e_ini(3) + e_ini(4)*e_ini(1)),(1-2*(e_ini(1)^2 + e_ini(2)^2)));

        B_body = eci2body(B_eci_ini, psi_ini, theta_ini, phi_ini);
        H_body = B_body/p.mu0;
        
        dH_dt_vector_ini = -cross(w0, H_body); % Calculate dH/dt at event
        
        
        if dH_dt_vector_ini(1) > 0 
            p.current_Hc_signs(1) = -1;
        else 
            p.current_Hc_signs(1) = 1;
        end

        if dH_dt_vector_ini(2) > 0 
            p.current_Hc_signs(2) = -1;
        else 
            p.current_Hc_signs(2) = 1;
        end
        
        if dH_dt_vector_ini(3) > 0 
            p.current_Hc_signs(3) = -1;
        else 
            p.current_Hc_signs(3) = 1;
        end


while t_current_segment_start < tspan(2)

    % Define ODE options for the current segment, including event function
   options = odeset('Events', @(t,state) sign_change_event(t, state, p), ...
                   'RelTol', 1e-7, 'AbsTol', 1e-7); % Or even tighter if needed

    % Options for higher precision and error control
    
    
    
    % Solve the ODE for the current segment
    [t_segment, state_segment, t_at_events, state_at_events, indices_of_rods_triggered] = ...
        ode45(@(t,state)physics(t,state,p), ... % Pass 'p' struct to motion function
              [t_current_segment_start, tspan(2)], ... % Integrate until end of full time span or event
              current_state_restart, options);

    % Append results from this segment to overall solution
    % Handle potential duplicates at segment boundaries for continuous plotting
    if isempty(all_t_sol) % First segment
        all_t_sol = t_segment;
        all_state_sol = state_segment;
        prev_last = last;
    else
        % Find where the new segment starts relative to previous accumulated data
        % This ensures no duplicate points if ode45 restarts exactly at previous end
        start_idx_in_segment = find(t_segment > all_t_sol(end), 1);
        if isempty(start_idx_in_segment) && ~isempty(t_segment) && t_segment(1) == all_t_sol(end)
            % If first point of new segment is identical to last point of old segment, skip it
            all_t_sol = [all_t_sol; t_segment(2:end)];
            all_state_sol = [all_state_sol; state_segment(2:end,:)];
        elseif ~isempty(start_idx_in_segment)
            % Normal case: append unique new points
            all_t_sol = [all_t_sol; t_segment(start_idx_in_segment:end)];
            all_state_sol = [all_state_sol; state_segment(start_idx_in_segment:end,:)];
        end
        prev_last = last;
    end

    last = length(unique(all_t_sol));
    for j = prev_last+1:last
        signs = [signs; p.current_Hc_signs];
    end

    prev_last = last;

    % Check if an event was detected
    % Check if an event was detected
    if ~isempty(t_at_events)
        % Update the start time and state for the next segment
        % Use the last event if multiple events happen exactly simultaneously
        t_current_segment_start = t_at_events(end);
        current_state_restart = state_at_events(end,:)'; % Transpose to match column vector state0
        
        % --- CRITICAL: Re-calculate dH/dt at the event time to determine correct Hc_signs ---
        % This ensures the sign is correct for the new segment, directly applying the rule.
        
        current_jd_event = p.jd + t_current_segment_start / (24 * 3600);
        gmst_rad_event = gmst_rad_from_jd(current_jd_event);
        
        % Extract necessary state variables from current_state_restart
        r_eci_event = current_state_restart(1:3);
        w_event = current_state_restart(11:13);

        % Recalculate B_body and H_body at the event point (necessary to get dH/dt)
        r_ecef = eci2ecef(r_eci_event, gmst_rad_event); 
        lat_long = ecef2lat_long(r_ecef);

        [Bn_nT, Be_nT, Bd_nT] = igrf("1-Jul-2025", lat_long(1)*(180/pi), lat_long(2)*(180/pi), (p.semi_major_orbit)/1000,'geocentric');
        Bu_nT = -Bd_nT;
        Bn = Bn_nT*1e-9;
        Be = Be_nT*1e-9;
        Bu = Bu_nT*1e-9;
        B_neu_ini = [Bn; Be; Bu];
        theta3 = 90 + atan2(r_ecef(2),r_ecef(1));
        theta1 = 90 - atan2(r_ecef(3),(sqrt(r_ecef(1)^2 + r_ecef(2)^2)));
        B_ecef_ini = neu2ecef(B_neu_ini,theta3,theta1);
        B_eci_ini = ecef2eci(B_ecef_ini, gmst_rad_event);

        e_ini = state_at_events(end,7:10);
        psi_ini = atan2(2*(e_ini(1)*e_ini(2) + e_ini(4)*e_ini(3)),(1-2*(e_ini(2)^2 + e_ini(3)^2)));
        theta_ini = asin(2*(e_ini(4)*e_ini(2) - e_ini(1)*e_ini(3)));
        phi_ini = atan2(2*(e_ini(2)*e_ini(3) + e_ini(4)*e_ini(1)),(1-2*(e_ini(1)^2 + e_ini(2)^2)));

        B_body = eci2body(B_eci_ini, psi_ini, theta_ini, phi_ini);
        H_body = B_body/p.mu0;
        
        dH_dt_vector_ini = -cross(w_event, H_body); % Calculate dH/dt at event
        
        
        if dH_dt_vector_ini(1) > 0 
            p.current_Hc_signs(1) = -1;
        else 
            p.current_Hc_signs(1) = 1;
        end

        if dH_dt_vector_ini(2) > 0 
            p.current_Hc_signs(2) = -1;
        else 
            p.current_Hc_signs(2) = 1;
        end
        
        if dH_dt_vector_ini(3) > 0 
            p.current_Hc_signs(3) = -1;
        else 
            p.current_Hc_signs(3) = 1;
        end


        fprintf('Time %.4f s: Hc_signs updated for restart: Y-rod: %.0f, Z-rod: %.0f (dH/dt_y=%.2e, dH/dt_z=%.2e)\n', ...
            t_current_segment_start, p.current_Hc_signs(1), p.current_Hc_signs(2), ...
            dH_dt_vector_ini(2), dH_dt_vector_ini(3));

    else
        % No event detected in this segment, means we've reached the end of tspan_full
        t_current_segment_start = tspan(2);
    end
end
fprintf('Simulation complete. Total time simulated: %.2f seconds.\n', all_t_sol(end));

[all_t_sol_unique, unique_indices] = unique(all_t_sol, 'stable');
all_state_sol_unique = all_state_sol(unique_indices, :);
% Extract translational state
r_eci_sol = all_state_sol_unique(:, 1:3);
v_eci_sol = all_state_sol_unique(:, 4:6);
% Extract rotational state
e_sol = all_state_sol_unique(:, 7:10);
angular_velocities_sol = all_state_sol_unique(:, 11:13);

% Convert quaternions to Euler angles (using atan2 for robustness)
% Your definition of euler parameters is [e1 e2 e3 eta]
psi_sol = atan2(2*(e_sol(:,1).*e_sol(:,2) + e_sol(:,4).*e_sol(:,3)), ...
               (1-2*(e_sol(:,2).^2 + e_sol(:,3).^2)));
theta_sol = asin(2*(e_sol(:,4).*e_sol(:,2) - e_sol(:,1).*e_sol(:,3)));
phi_sol = atan2(2*(e_sol(:,2).*e_sol(:,3) + e_sol(:,4).*e_sol(:,1)), ...
               (1-2*(e_sol(:,1).^2 + e_sol(:,2).^2)));


% Recalculate Bx_I,By_I,Bz_I, Bx_b, By_b, Bz_b for plotting purposes (only once after sim is complete)
Bx_I = zeros(length(all_t_sol_unique), 1);
By_I = zeros(length(all_t_sol_unique), 1);
Bz_I = zeros(length(all_t_sol_unique), 1);
Bx_b = zeros(length(all_t_sol_unique), 1);
By_b = zeros(length(all_t_sol_unique), 1);
Bz_b = zeros(length(all_t_sol_unique), 1);
Bhyst_x = zeros(length(all_t_sol_unique), 1);
Bhyst_y = zeros(length(all_t_sol_unique), 1);
Bhyst_z = zeros(length(all_t_sol_unique), 1);
Hx_b = zeros(length(all_t_sol_unique), 1);
Hy_b = zeros(length(all_t_sol_unique), 1);
Hz_b = zeros(length(all_t_sol_unique), 1);

for i = 1:length(all_t_sol_unique)
    current_t_sim = all_t_sol_unique(i);
    e = e_sol(i,:)'; % Quaternion for current time step [e1; e2; e3; eta]
    current_r_eci = r_eci_sol(i,:)';     % Position for current time step

    % Get current Julian Date and GMST
    current_jd = p.jd + current_t_sim / (24 * 3600);
    gmst_rad = gmst_rad_from_jd(current_jd);

    % Convert ECI position to ECEF and then to LLA
    r_ecef = eci2ecef(current_r_eci, gmst_rad);
    lat_long = ecef2lat_long(r_ecef);

    % Corrected IGRF altitude input: use current norm(r_eci)
    [Bn_nT, Be_nT, Bd_nT] = igrf("1-Jul-2025", lat_long(1)*(180/pi), lat_long(2)*(180/pi), (p.semi_major_orbit)/1000,'geocentric');
    Bu_nT = -Bd_nT;
    Bn = Bn_nT*1e-9;
    Be = Be_nT*1e-9;
    Bu = Bu_nT*1e-9;
    B_neu = [Bn; Be; Bu];
    theta3 = 90 + atan2(r_ecef(2),r_ecef(1));
    theta1 = 90 - atan2(r_ecef(3),sqrt(r_ecef(1)^2 + r_ecef(2)^2));
    B_ecef = neu2ecef(B_neu,theta3,theta1);
    B_eci = ecef2eci(B_ecef, gmst_rad);

    psi = atan2(2*(e(1)*e(2) + e(4)*e(3)),(1-2*(e(2)^2 + e(3)^2)));
    theta = asin(2*(e(4)*e(2) - e(1)*e(3)));
    phi = atan2(2*(e(2)*e(3) + e(4)*e(1)),(1-2*(e(1)^2 + e(2)^2)));

    B_body = eci2body(B_eci, psi, theta, phi);
    H_body = B_body/p.mu0;
    B_hyst = (2/pi)*p.Bs*atan(p.p0*(H_body + p.Hc*signs(i,:)'));

    Bhyst_x(i) = B_hyst(1);
    Bhyst_y(i) = B_hyst(2);
    Bhyst_z(i) = B_hyst(3);
    Hx_b(i) = H_body(1);
    Hy_b(i) = H_body(2);
    Hz_b(i) = H_body(3);
  
end

% Plotting
figure;
plot3(r_eci_sol(:,1),r_eci_sol(:,2),r_eci_sol(:,3), "b", LineWidth = 2)
grid on
hold on
[X ,Y ,Z] = sphere;
X = X*p.radius_earth;
Y = Y*p.radius_earth;
Z = Z*p.radius_earth;
surf(X,Y,Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none') % Make Earth transparent
axis equal;
title('Satellite Orbit');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

figure;
plot(all_t_sol_unique, Bx_I * 1e9, "b-", LineWidth = 2); % Convert to nT
grid on
hold on
plot(all_t_sol_unique, By_I * 1e9, "g-", LineWidth = 2);
plot(all_t_sol_unique, Bz_I * 1e9, "r-", LineWidth = 2);
xlabel("Time (sec)");
ylabel("B magnitude (nT)");
legend("B_x ECI","B_y ECI","B_z ECI", 'Location', 'best');
title('Magnetic Field in ECI Frame');

figure;
plot(all_t_sol_unique, Bx_b * 1e9, "b-", LineWidth = 2); % Convert to nT
grid on
hold on
plot(all_t_sol_unique, By_b * 1e9, "g-", LineWidth = 2);
plot(all_t_sol_unique, Bz_b * 1e9, "r-", LineWidth = 2);
xlabel("Time (sec)");
ylabel("B magnitude (T)");
legend("B_x Body","B_y Body","B_z Body", 'Location', 'best');
title('Magnetic Field in Body Frame');

figure;
subplot(3,2,1);
plot(all_t_sol_unique,psi_sol * 180/pi)
grid on
title("Psi (Yaw) vs Time (deg)")
xlabel("Time (sec)")
ylabel("Psi (deg)")
subplot(3,2,3);
plot(all_t_sol_unique,theta_sol * 180/pi)
grid on
title("Theta (Pitch) vs Time (deg)")
xlabel("Time (sec)")
ylabel("Theta (deg)")
subplot(3,2,5);
plot(all_t_sol_unique,phi_sol * 180/pi)
grid on
title("Phi (Roll) vs Time (deg)")
xlabel("Time (sec)")
ylabel("Phi (deg)")
subplot(3,2,2);
plot(all_t_sol_unique,angular_velocities_sol(:,1) * 180/pi)
grid on
title("Angular Velocity \omega_x vs Time (deg/s)")
xlabel("Time (sec)")
ylabel("\omega_x (deg/s)")
subplot(3,2,4);
plot(all_t_sol_unique,angular_velocities_sol(:,2) * 180/pi)
grid on
title("Angular Velocity \omega_y vs Time (deg/s)")
xlabel("Time (sec)")
ylabel("\omega_y (deg/s)")
subplot(3,2,6);
plot(all_t_sol_unique,angular_velocities_sol(:,3) * 180/pi)
grid on
title("Angular Velocity \omega_z vs Time (deg/s)")
xlabel("Time (sec)")
ylabel("\omega_z (deg/s)")
sgtitle('Satellite Attitude Dynamics & Rotation');

figure;
subplot(3,1,1);
plot(Hx_b,Bhyst_x);
grid on
title("Bhyst_x vs Hx_b")
xlabel("Hx_b")
ylabel("Bhyst_x")

subplot(3,1,2);
plot(Hy_b,Bhyst_y);
grid on
title("Bhyst_y vs Hy_b")
xlabel("Hy_b")
ylabel("Bhyst_y")

subplot(3,1,3);
plot(Hz_b,Bhyst_z);
grid on
title("Bhyst_z vs Hz_b")
xlabel("Hz_b")
ylabel("Bhyst_z")