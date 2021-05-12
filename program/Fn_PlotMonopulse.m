function Fn_PlotMonopulse(theta_s,theta3,az_beam,el_beam)

a=theta3;
plot(az_beam,el_beam,'k+','LineWidth',2);
axis([az_beam-a az_beam+a el_beam-a el_beam+a]), 
axis square, hold on;
plot([az_beam-a az_beam+a],[el_beam,el_beam],'Color',[0.7 0.7 0.7]);
plot([az_beam az_beam],[el_beam-a,el_beam+a],'Color',[0.7 0.7 0.7]);
xlabel('azimuth (deg)'), ylabel('elevation (deg)')
title('Angle Scope');
% text(az_beam-theta_s,el_beam-theta3+theta3/5,'Monopulse','Color','black','FontSize',14); 

bm=theta_s/sqrt(2); % beam_max
plot(az_beam+bm,el_beam+bm,'ko','MarkerFace','k','MarkerSize',3)
cx1=bm-theta3/2:0.01:bm+theta3/2;
plot(az_beam+cx1,el_beam+bm+sqrt((theta3/2)^2-(cx1-bm).^2),'Color',[0.7 0.7 0.7]), 
plot(az_beam+cx1,el_beam+bm-sqrt((theta3/2)^2-(cx1-bm).^2),'Color',[0.7 0.7 0.7]),
plot(az_beam-bm,el_beam+bm,'ko','MarkerFace','k','MarkerSize',3)
cx2=-bm-theta3/2:0.01:-bm+theta3/2;
plot(az_beam+cx2,el_beam+bm+sqrt((theta3/2)^2-(cx2+bm).^2),'Color',[0.7 0.7 0.7]'), 
plot(az_beam+cx2,el_beam+bm-sqrt((theta3/2)^2-(cx2+bm).^2),'Color',[0.7 0.7 0.7]),
plot(az_beam-bm,el_beam-bm,'ko','MarkerFace','k','MarkerSize',3)
cx3=-bm-theta3/2:0.01:-bm+theta3/2;
plot(az_beam+cx3,el_beam-bm+sqrt((theta3/2)^2-(cx3+bm).^2),'Color',[0.7 0.7 0.7]), 
plot(az_beam+cx3,el_beam-bm-sqrt((theta3/2)^2-(cx3+bm).^2),'Color',[0.7 0.7 0.7]),
plot(az_beam+bm,el_beam-bm,'ko','MarkerFace','k','MarkerSize',3)
cx4=bm-theta3/2:0.01:bm+theta3/2;
plot(az_beam+cx4,el_beam-bm+sqrt((theta3/2)^2-(cx4-bm).^2),'Color',[0.7 0.7 0.7]), 
plot(az_beam+cx4,el_beam-bm-sqrt((theta3/2)^2-(cx4-bm).^2),'Color',[0.7 0.7 0.7]), 
