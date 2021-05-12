function Fn_PlotSearchScope(max_range,target,BT_range)
% max_range [km]
% beam_width [deg]
% target = (x_target, y_target) [km]
x_target=target(1); y_target=target(2);
cx=-max_range:0.01:max_range;
plot(cx,sqrt(max_range^2-cx.^2),'Color',[0.7 0.7 0.7]); 
axis([-(max_range+20) (max_range+20) -(max_range+20) (max_range+20)]), 
axis square, hold on;
plot(cx,-sqrt(max_range^2-cx.^2),'Color',[0.7 0.7 0.7]);
radius=max_range/3;
cx=-radius:0.01:radius;
plot(cx,sqrt(radius^2-cx.^2),'Color',[0.7 0.7 0.7]);
plot(cx,-sqrt(radius^2-cx.^2),'Color',[0.7 0.7 0.7]);
radius=max_range/3*2;
cx=-radius:0.01:radius;
plot(cx,sqrt(radius^2-cx.^2),'Color',[0.7 0.7 0.7]);
plot(cx,-sqrt(radius^2-cx.^2),'Color',[0.7 0.7 0.7]);
radius=BT_range;
cx=-radius:0.01:radius;
plot(cx,sqrt(radius^2-cx.^2),'Color','red');
plot(cx,-sqrt(radius^2-cx.^2),'Color','red');
plot([-max_range-20 max_range+20],[0 0],'Color',[0.7 0.7 0.7]);
plot([0 0],[-max_range-20 max_range+20],'Color',[0.7 0.7 0.7]);

plot([max_range max_range],[0 -(max_range+20)],':','Color',[0.3 0.3 0.3]);
instr=[num2str(max_range),'km'];
str={'Maximum','Range',instr};
text(max_range+2,-max_range+8,str,'Color','black');
plot(0,0,'k+','LineWidth',1);
plot(x_target,y_target,'rx','LineWidth',2);
plot([x_target -max_range],[y_target max_range],':','Color','red')
text(-max_range-10,max_range+5,'target','Color','red');
title('Radar Scope');
xlabel('range (km)');
ylabel('range (km)');