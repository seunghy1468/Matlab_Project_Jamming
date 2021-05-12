function result=Fn_TargetMove(ip,vel,ad,num,intv)
% ip : initial position (x,y) [km]
% h: height of target (fixed) [km]
% vel : velocity of target [m/s]
% ad : angle of direction [deg]
% num : number of time-line
% intv : time interval
x=ip(1); y=ip(2);
result=[x,y];
for k=1:num-1,
    x=x+intv*vel*cos(ad*pi/180)/1000;
    y=y+intv*vel*sin(ad*pi/180)/1000;
    result=[result; x,y];
end
