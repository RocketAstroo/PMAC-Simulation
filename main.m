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
  wx0 = 0.6;
  wy0 = 0.6;
  wz0 = 0.6;
  w0 = [wx0; wy0; wz0];

%Initialization of date, time 
start_jd = gregorian_to_julian(p.date); % Julian Date of simulation start

%Step size of the RK4 step fixed solver and timespan of the simulation
h = 1;
no_of_orbits = 175;
tspan = [0, 5000];
time_stamps = 0:h:tspan(2);
len = length(time_stamps);
current_pt = 1;

 psi_orb = zeros(len, 1);
  theta_orb = zeros(len, 1);
  phi_orb = zeros(len, 1);

  psi_eci = zeros(len, 1);
  theta_eci = zeros(len, 1);
  phi_eci = zeros(len, 1);


%Initialization of the overall state vector and dhdt values
  state = zeros(length(r0_eci)+length(v0_eci)+length(e0)+length(w0),len);
  state(:,1) = [r0_eci;v0_eci;e0;w0];

%Initialize Hc signs
current_jd = p.jd + tspan(1)/(24*3600);
gmst_rad = gmst_rad_from_jd(current_jd);
    
r_ecef = eci2ecef(state(1:3,1), gmst_rad); 
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

        e_ini = state(7:10,1);
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

  
  dHdt = zeros(3,len);
  dHdt(:,1) =  dH_dt_vector_ini;

    for j = current_pt:len-1

        k1 = physics(time_stamps(j),state(:,j),p);
        k2 = physics(time_stamps(j) + (h/2),state(:,j) + (h/2)*k1,p);
        k3 = physics(time_stamps(j) + (h/2),state(:,j) + (h/2)*k2,p);
        k4 = physics(time_stamps(j) + h,state(:,j) + h*k3,p);
        
        state(:,j+1) = state(:,j) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

        e1 = state(7,j+1);
        e2 = state(8,j+1);
        e3 = state(9,j+1);
        e4 = state(10,j+1);

     psi_eci(j+1) = atan(2*(e1*e2 + e4*e3)/(1-2*(e2^2 + e3^2)));
    theta_eci(j+1) = asin(2*(e4*e2 - e1*e3));
      phi_eci(j+1) = atan(2*(e2*e3 + e4*e1)/(1-2*(e1^2 + e2^2)));

       angle_wrt_orb = eci2orb(state(1:3,j+1),state(4:6,j+1), psi_eci(j+1), theta_eci(j+1), phi_eci(j+1));
       psi_orb(j+1) = angle_wrt_orb(1);
       theta_orb(j+1) = angle_wrt_orb(2);
       phi_orb(j+1) = angle_wrt_orb(3);

        dHdt(:, j+1) = dHdt_val(time_stamps(j+1),state(:,j+1),p);

        if dHdt(1,j+1) > 0 
            p.current_Hc_signs(1) = -1;
        else 
            p.current_Hc_signs(1) = 1;
        end

        if dHdt(2,j+1) > 0 
            p.current_Hc_signs(2) = -1;
        else 
            p.current_Hc_signs(2) = 1;
        end
        
        if dHdt(3,j+1) > 0 
            p.current_Hc_signs(3) = -1;
        else 
            p.current_Hc_signs(3) = 1;
        end

    end

    disp("Simulation ended");


% Extract translational state
r_eci_sol = state(1:3,:);
v_eci_sol = state(4:6,:);
% Extract rotational state
e_sol = state(7:10, :);
angular_velocities_sol = state(11:13, :);

% Convert quaternions to Euler angles (using atan2 for robustness)
% Your definition of euler parameters is [e1 e2 e3 eta]
psi_sol = atan2(2*(e_sol(1, :).*e_sol(2, :) + e_sol(4, :).*e_sol(3, :)), ...
               (1-2*(e_sol(2, :).^2 + e_sol(3, :).^2)));
theta_sol = asin(2*(e_sol(4, :).*e_sol(2, :) - e_sol(1, :).*e_sol(3, :)));
phi_sol = atan2(2*(e_sol(2, :).*e_sol(3, :) + e_sol(4, :).*e_sol(1, :)), ...
               (1-2*(e_sol(1, :).^2 + e_sol(2, :).^2)));


% Recalculate Bx_I,By_I,Bz_I, Bx_b, By_b, Bz_b for plotting purposes (only once after sim is complete)
Bx_I = zeros(len, 1);
By_I = zeros(len, 1);
Bz_I = zeros(len, 1);
Bx_b = zeros(len, 1);
By_b = zeros(len, 1);
Bz_b = zeros(len, 1);
Bhyst_x = zeros(len, 1);
Bhyst_y = zeros(len, 1);
Bhyst_z = zeros(len, 1);
Hx_b = zeros(len, 1);
Hy_b = zeros(len, 1);
Hz_b = zeros(len, 1);

for i = 1:len
    current_t_sim = time_stamps(i);
    e = e_sol(:,i); % Quaternion for current time step [e1; e2; e3; eta]
    current_r_eci = r_eci_sol(:,i);     % Position for current time step

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
    Bx_I(i) = B_eci(1);
    By_I(i) = B_eci(2);
    Bz_I(i) = B_eci(3);

    psi = atan2(2*(e(1)*e(2) + e(4)*e(3)),(1-2*(e(2)^2 + e(3)^2)));
    theta = asin(2*(e(4)*e(2) - e(1)*e(3)));
    phi = atan2(2*(e(2)*e(3) + e(4)*e(1)),(1-2*(e(1)^2 + e(2)^2)));

    B_body = eci2body(B_eci, psi, theta, phi);
    Bx_b(i) = B_body(1);
    By_b(i) = B_body(2);
    Bz_b(i) = B_body(3);

    H_body = B_body/p.mu0;
    
    dH_dt = -cross(state(11:13,i), H_body);

        if dH_dt(1) > 0 
            p.current_Hc_signs(1) = -1;
        else 
            p.current_Hc_signs(1) = 1;
        end

        if dH_dt(2) > 0 
            p.current_Hc_signs(2) = -1;
        else 
            p.current_Hc_signs(2) = 1;
        end
        
        if dH_dt(3) > 0 
            p.current_Hc_signs(3) = -1;
        else 
            p.current_Hc_signs(3) = 1;
        end


    B_hyst = (2/pi)*p.Bs*atan(p.p0*(H_body + p.Hc*p.current_Hc_signs));

    Bhyst_x(i) = B_hyst(1);
    Bhyst_y(i) = B_hyst(2);
    Bhyst_z(i) = B_hyst(3);
    Hx_b(i) = H_body(1);
    Hy_b(i) = H_body(2);
    Hz_b(i) = H_body(3);
  
end

% Plotting
figure;
plot3(r_eci_sol(1,:),r_eci_sol(2, :),r_eci_sol(3, :), "b", LineWidth = 2)
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
plot(time_stamps, Bx_I * 1e9, "b-", LineWidth = 2); % Convert to nT
grid on
hold on
plot(time_stamps, By_I * 1e9, "g-", LineWidth = 2);
plot(time_stamps, Bz_I * 1e9, "r-", LineWidth = 2);
xlabel("Time (sec)");
ylabel("B magnitude (nT)");
legend("B_x ECI","B_y ECI","B_z ECI", 'Location', 'best');
title('Magnetic Field in ECI Frame');

figure;
plot(time_stamps, Bx_b * 1e9, "b-", LineWidth = 2); % Convert to nT
grid on
hold on
plot(time_stamps, By_b * 1e9, "g-", LineWidth = 2);
plot(time_stamps, Bz_b * 1e9, "r-", LineWidth = 2);
xlabel("Time (sec)");
ylabel("B magnitude (T)");
legend("B_x Body","B_y Body","B_z Body", 'Location', 'best');
title('Magnetic Field in Body Frame');

figure;
subplot(3,2,1);
plot(time_stamps, psi_orb * 180/pi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Actual $\psi$');
hold on
grid on;
title("Yaw ($\psi$) vs Time", 'Interpreter', 'latex');
xlabel("Time (sec)");
ylabel("Angle (deg)");
hold off;

subplot(3,2,3);
plot(time_stamps, theta_orb * 180/pi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Actual $\theta$');
hold on
grid on;
title("Pitch ($\theta$) vs Time", 'Interpreter', 'latex');
xlabel("Time (sec)");
ylabel("Angle (deg)");
hold off;

subplot(3,2,5);
plot(time_stamps, phi_orb * 180/pi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Actual $\phi$');
hold on

grid on;
title("Roll ($\phi$) vs Time", 'Interpreter', 'latex');
xlabel("Time (sec)");
ylabel("Angle (deg)");
hold off;

subplot(3,2,2);
plot(time_stamps,angular_velocities_sol(1, :) * 180/pi)
grid on
title("Angular Velocity \omega_x vs Time (deg/s)")
xlabel("Time (sec)")
ylabel("\omega_x (deg/s)")
subplot(3,2,4);
plot(time_stamps,angular_velocities_sol(2, :) * 180/pi)
grid on
title("Angular Velocity \omega_y vs Time (deg/s)")
xlabel("Time (sec)")
ylabel("\omega_y (deg/s)")
subplot(3,2,6);
plot(time_stamps,angular_velocities_sol(3, :) * 180/pi)
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
