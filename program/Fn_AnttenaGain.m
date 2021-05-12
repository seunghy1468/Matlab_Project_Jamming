function G=Fn_AnttenaGain(g0, theta3, theta_s, azi_angle, elv_angle)
bm=theta_s/sqrt(2); % beam_max

ang_dist=sqrt((bm-azi_angle)^2+(bm-elv_angle)^2);
g1=g0*exp(-2*log(2)*(ang_dist/theta3).^2);

ang_dist=sqrt((-bm-azi_angle)^2+(bm-elv_angle)^2);
g2=g0*exp(-2*log(2)*(ang_dist/theta3).^2);

ang_dist=sqrt((-bm-azi_angle)^2+(-bm-elv_angle)^2);
g3=g0*exp(-2*log(2)*(ang_dist/theta3).^2);

ang_dist=sqrt((bm-azi_angle)^2+(-bm-elv_angle)^2);
g4=g0*exp(-2*log(2)*(ang_dist/theta3).^2);

G=[g1 g2 g3 g4];