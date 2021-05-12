function [fd, direction]=Fn_DopplerFreq(target,h_radar,h_target,ad_target,vel,f_RF)
% f_RF [MHz]
% vel [m/s]
c=3; % *10^8 [m/s]
direction=(-target(1)*cos(ad_target*pi/180)-target(2)*sin(ad_target*pi/180))...
    /sqrt(target(1)^2+target(2)^2+(h_radar-h_target)^2);
fd=2*direction*vel*f_RF/c*10^(-2); % [Hz]
% result=[fd, direction];