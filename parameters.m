%This file contains a function called parameters which returns a structure
%p which contains all the necessary parameter values (all are in SI units
%and radians until unless specially mentioned)

function p = parameters()

%Constants 
    p.mu0 = 4 * pi * 10^(-7);
    p.gravitational_const = 6.6743e-11;

% Earth parameters
    p.mass_earth = 5.9722e24; 
    p.radius_earth = 6.3781e6; 
    p.w_earth = 7.292115e-5; 

% Satellite parameters
    p.mass_sat = 3;
    p.square_side_sat = 10e-2; 
    p.length_sat = 34.05e-2; 
    p.Ixx = (p.mass_sat/12) * (p.square_side_sat^2 + p.square_side_sat^2); %X principal axis is the long axis of the sat(prolate body)
    p.Iyy = (p.mass_sat/12) * (p.length_sat^2 + p.square_side_sat^2);
    p.Izz = p.Iyy;
    p.Inertia_tensor_sat = diag([p.Ixx,p.Iyy,p.Izz]);

%Bar magnet parameters (kept along X-axis of the body frame of the sat)
    p.moment_magnet = 0.3; 
    p.moment_magnet_vec = [p.moment_magnet; 0; 0];

%Hysterisis rod design parameters (kept along Y and Z axes of the body frame of the sat)
    p.Rod_length = 0.095; 
    p.Rod_diameter = 0.001; 
    p.Hc = 12; 
    p.Br = 0.004; 
    p.Bs = 0.027;
    p.rod_vol = p.Rod_length * (pi * (p.Rod_diameter / 2)^2); 
    p.no_rods = [3;3;3]; %along X,Y and Z 
    p.p0 = (1/p.Hc)*tan((pi*p.Br)/(2*p.Bs));
    p.current_Hc_signs = [-1;-1;-1]; %sign of Hc in Flatley & Henretty model of the hysterisis loop

%Orbital parameters
    p.inclination_orbit = 98 * (pi / 180); %SSPO orbit
    p.altitude_orbit = 600e3; 
    p.semi_major_orbit = p.altitude_orbit + p.radius_earth;
    p.vel_orbit = sqrt((p.gravitational_const * p.mass_earth) / (p.semi_major_orbit)); 
    p.period_orbit = 2 * pi * (p.semi_major_orbit) / p.vel_orbit; 

%Date
    p.date = [2025, 07, 2, 4, 47, 0]; % Y, M, D, H, M, S (IST) - Use current time for consistency
    p.jd = gregorian_to_julian(p.date); % Julian Date of simulation start

end