function rho_half=Fn_PulseSpread(elev_rad,psi_g,R,G,theta3)
c=3;
max_R=R/sin(pi-2*psi_g)*(sin(psi_g+elev_rad)+sin(psi_g-elev_rad)); % max distance of reflection, km
rho_max=(max_R-R)*10/c; % maximum delay spread, micro sec (*10^{-6})
ref_rate=sqrt((G(1)*sin(psi_g))^2+(G(2)*cos(psi_g))^2); % reflection amplitude
antenna_gain=exp(-(elev_rad+psi_g)^2/((theta3*pi/180)^2/log(2))); % antenna_gain
rho_half=(1-sqrt(0.5))*rho_max/(1-antenna_gain*ref_rate); % 3dB delay spread, micro sec
