function result=Fn_TargetPoSearch(target,h,beam_width,R_est)
% for search radar
x_target=target(1);
y_target=target(2);
r_target=sqrt(x_target^2+y_target^2+h^2);
if x_target>0 && y_target==0, a_target=0;
elseif x_target>0 && y_target>0, a_target=atan(y_target/x_target);
elseif x_target==0 && y_target>0, a_target=pi/2;
elseif x_target<0 && y_target>0, a_target=pi+atan(+y_target/x_target);
elseif x_target<0 && y_target==0, a_target=pi;
elseif x_target<0 && y_target<0, a_target=pi+atan(y_target/x_target);
elseif x_target==0 && y_target<0, a_target=3*pi/2;
else a_target=2*pi+atan(y_target/x_target);   
end

angle_bin=360/beam_width; % [deg]
for k=0:angle_bin-1,
    if a_target<=beam_width*pi/180*(k+1), 
        a_est=beam_width*pi/180*(2*k+1)/2; break;
    end
end
x_est=R_est*cos(a_est);
y_est=R_est*sin(a_est);
result=[x_est y_est];