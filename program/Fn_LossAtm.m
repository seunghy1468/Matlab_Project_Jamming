function gamma=Fn_LossAtm(freq,height)
% Van Vleck Equation
% freq : 레이더 주파수, GHz
% height : 수신기 높이, km
v1=0.018; v2=0.05;
lambda=3e+10./(freq*1e+9); % wavelength, cm
T=288-(6.7*height); % 15deg, temperature at different heights, Kelvin deg
P=1015*((1-0.02257*height).^(5.2561)); % air pressure at different heights
%---- gamma_O2 % dB/km
A1=(0.4909*(P.^2)*v1)./(T.^(5/2));
A2=1+(2.904e-4*(lambda^2)*(P.^2)*(v1^2)./T);
A3=1+(0.5*v2)/((lambda^2)*v1);
gamma_O2=(A1./A2)*A3;
%---- gamma_H2O % dB/km
for k=1:length(height),
if height(k)<0.381, wvd(k)=6.18; % wvd = water vapor density, g/m^3
elseif height(k)<1.143, wvd(k)=4.93;
elseif height(k)<2.286, wvd(k)=3.74;
elseif height(k)<4.572, wvd(k)=2.01;
elseif height(k)<7.620, wvd(k)=0.34;
elseif height(k)<10.668, wvd(k)=0.05;
else wvd(k)=0;
end
end
B1=1.852*(3.165e-6)*wvd.*(P.^2)./(T.^(3/2));
B2=3.43/(lambda^2);
B3=(1-0.742*lambda)^2;
B4=(1+1.742*lambda)^2;
B5=(2.853e-6)*(lambda^2)*P.^2./T;
gamma_H2O=B1.*(1./(B3+B5)+1./(B4+B5)+B2);
%--- total gamma
gamma=gamma_O2+gamma_H2O;