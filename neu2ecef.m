function v_ecef = neu2ecef(v_neu,theta3,theta1) 

C = [-cos(theta1)*sin(theta3), sin(theta1)*sin(theta3), cos(theta3); ...
    cos(theta1)*cos(theta3), -sin(theta1)*cos(theta3), sin(theta3); ...
    sin(theta1), cos(theta1), 0];

v_ecef = C*v_neu;

end