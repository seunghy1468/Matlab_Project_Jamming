azi_angle=2; % target azi angle from beam center, deg
elv_angle=-2; % target elv angle from beam center, deg
theta_s=5; % deg
theta3=10; % deg
g0=1;
km=2*sqrt(2)*log(2)*theta_s/theta3; % antenna constant
V=Fn_AnttenaGain(g0, theta3, theta_s, azi_angle, elv_angle);

%%--- Monopulse Diff and Sum (Angle Estimation)
% target azimuth and elevation: azim, elev [deg]
% beam center: az_beam, el_beam [deg]
% target angles from beam center: azi_angle, elv_angle [deg]
% estimated angles from beam center: hat_azi, hat_elv
delta_az=(V(1)-V(2))-(V(3)-V(4));
delta_el=(V(1)+V(2))-(V(3)+V(4));
hat_azi=delta_az*theta3/(4*km) % estimation of azi_angle
hat_elv=delta_el*theta3/(4*km) % estimation of elv_angle
    