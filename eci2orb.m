function angle_wrt_orb = eci2orb(r,v, psi, theta, phi)

z_orb = r/norm(r);
y_orb = cross(r,v)/norm(cross(r,v));
x_orb = cross(y_orb,z_orb);

C_eci2orb = [x_orb';y_orb';z_orb'];

C_eci2body = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta); ...
    (-cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi)), (cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi)), sin(phi)*cos(theta); ...
    (sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi)), (-sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)), cos(phi)*cos(theta)];

%C_eci2body = C_orb2body*C_eci2orb
%C_eci2body*C_eci2orb' = C_orb2body

C_orb2body = C_eci2body*C_eci2orb';

%ψ=atan2(C12​,C11​)
psi_orb = atan2(C_orb2body(1,2),C_orb2body(1,1));

%θ=arcsin(−C13​)
theta_orb = asin(-C_orb2body(1,3));

%ϕ=atan2(C23​,C33​)
phi_orb = atan2(C_orb2body(2,3), C_orb2body(3,3));

angle_wrt_orb = [psi_orb;theta_orb; phi_orb];




end