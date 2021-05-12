function Fn_PlotAngScope(theta3,est_bin_az,est_bin_el,signal_miss,flag)
xc=mod(est_bin_az+theta3/2,360);
yc=est_bin_el+theta3/2;

plot(xc,yc,'k+','LineWidth',2) %,'MarkerFace','k','MarkerSize',3)
a=theta3;
axis([xc-a xc+a yc-a yc+a]), axis square, hold on;
xlabel('azimuth (deg)'), ylabel('elevation (deg)')
title('Angle Scope');
plot(xc+[-a a],yc+[0 0],'Color',[0.7 0.7 0.7]);
plot(xc+[0 0],yc+[-a a],'Color',[0.7 0.7 0.7]);
cx1=-theta3/2:0.01:+theta3/2;
    
if signal_miss==0 && flag==1
    plot(xc+cx1,yc+sqrt((theta3/2)^2-(cx1).^2),'Color','blue')
    plot(xc+cx1,yc-sqrt((theta3/2)^2-(cx1).^2),'Color','blue')
else    
    plot(xc+cx1,yc+sqrt((theta3/2)^2-(cx1).^2),'Color',[0.7 0.7 0.7])
    plot(xc+cx1,yc-sqrt((theta3/2)^2-(cx1).^2),'Color',[0.7 0.7 0.7])
end