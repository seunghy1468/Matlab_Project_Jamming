function [R,azim,elev]=Fn_RangeAzEl(target,h_target)

R=sqrt(target(1)^2+target(2)^2+h_target^2); % [km] target range
elev_rad=asin(h_target/R); % [rad] elevation angle
elev=elev_rad*180/pi; % elevation angle, rad
if target(1)>0 && target(2)>0, 
    azim=asin(target(2)/(R*cos(elev_rad)))*180/pi; % [deg] azimuth angle
elseif target(1)<0 && target(2)>0, 
    azim=180-asin(target(2)/(R*cos(elev_rad)))*180/pi; % [deg] azimuth angle
elseif target(1)<0 && target(2)<0,
    azim=180-asin(target(2)/(R*cos(elev_rad)))*180/pi; % [deg] azimuth angle
else  
    azim=360+asin(target(2)/(R*cos(elev_rad)))*180/pi; % [deg] azimuth angle    
end
% result=[R,azim,elev];