function result=Fn_PulseAmp(Ptotal,PRI,PW,R,f_RF,RCS,h_target,rainfall)
% f_RF [MHz]
c=3;  
% lambda=c*10^8/(f_RF*10^6); % [m] wavelength of radar signal
% A_factor= 10*log10(lambda^2*RCS / (4*pi)^3 / (R*1000)^4); % [dB] freq, RCS, range
Lfs=10*log10((4*pi)^3*(R*1000)^4*(f_RF*10^6/c*10^-8)^2/RCS); % [dB] two-way freespace loss
Latm=2*R*Fn_LossAtm(f_RF/1000, h_target); % [dB] atmospheric (O2,H2O) loss
Lrain=2*R*0.0002*(f_RF/1000)^2.25*rainfall; % [dB] rain loss
PP_rcv=Ptotal-(Lfs+Latm+Lrain); % [dB-W] received power
% PP_rcv=Ptotal-(Latm+Lrain)+A_factor; % [dB-W] received power
% PP_rcv=Ptotal+A_factor; % [dB-W] received power
result=sqrt(10^(PP_rcv/10)*2); % [V]

