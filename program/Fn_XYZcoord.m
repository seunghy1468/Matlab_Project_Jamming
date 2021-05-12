function result=Fn_XYZcoord(R,azim,elev,h_radar)
az=azim*pi/180;
el=elev*pi/180;
x=R*cos(el)*cos(az);
y=R*cos(el)*sin(az);
z=R*sin(el)+h_radar;
result=[x,y,z]; % [km]