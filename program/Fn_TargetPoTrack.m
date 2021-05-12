function result=Fn_TargetPoTrack(R,azi,elv)
% for tracking radar
x=R*cos(elv*pi/180)*cos(azi*pi/180);
y=R*cos(elv*pi/180)*sin(azi*pi/180);
z=R*sin(elv*pi/180);
result=[x y z];