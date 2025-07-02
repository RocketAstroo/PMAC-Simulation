function lat_long = ecef2lat_long(r_ecef)

lat = atan(r_ecef(3)/(sqrt(r_ecef(1)^2 + r_ecef(2)^2)));
long = atan(r_ecef(2)/r_ecef(1));
lat_long = [lat;long];

end