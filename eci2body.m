function v_body = eci2body(v_eci, psi, theta, phi)

C = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta); ...
    (-cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi)), (cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi)), sin(phi)*cos(theta); ...
    (sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi)), (-sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)), cos(phi)*cos(theta)];

v_body = C*v_eci;

end