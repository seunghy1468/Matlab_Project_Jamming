function Fn_RadarSearchContWave()

load each_sim.mat
ccc=1;

%% ===== 변수 계산 및 세팅 =====

%%--- time and freq
max_time=PRI; % [micro seq] duration per a simulation 
tt=1/Fs; % [micro sec] sampling interval in time
t=0:tt:max_time-tt; % [micro sec] horizontal axis in time
Nt=length(t); 
txwave=0:tt:PW-tt; % [micro sec] pulse (waveform) time
Ntx=length(txwave);
% ff=Fs/Nt; % [MHz/sample] sampling interval in freq
% f=(0:Nt-1)*ff; % [MHz] horizontal axis in freq

%%--- power, amplitude, threshold, and noise
PA = sqrt(2*10^((PP-30)/10)); % [V] pulse amplitude, sinusoidal waveform
SNR = ( sqrt(-log(Pfa)) - erfcinv(2*Pd) )^2 - 0.5;
lambda = c*10^8 / (f_RF*10^6); % [m] wavelength of radar signal
max_range = PRI*c/2/10; % [km] maximum range of radar
Pr_min = (10^((PP-30+Gt+Gr)/10)*(lambda^2)*RCS) / ...
    ((4*pi)^3*((max_range*1.05)*10^3)^4); % min rx power (calculated using the max range)
sigma_square = Pr_min/SNR; % noise variance
sigma = sqrt(sigma_square)/tt;
Vt = sqrt(2 * sigma_square * log(1/Pfa)); % matched filter threshold
% %--- for Combined Simulator (5F Lab)
% PA=5000; % [V] pulse amplitude
% PP=0.5*PW*PA^2/PRI/1000; % [kW] pulse power

%%--- IF
IF_result=Fn_IFcoeff(wav,amp_op,f_IF/1000,IF_bw,IF_type,IF_perf,PW,Fs);
la=IF_result(1); lb=IF_result(1+la+1);
a=IF_result(2:la+1);
b=IF_result(la+3:la+lb+2);
[transfer, w]=freqz(b, a, Ntx, Fs); 
ir=ifft(transfer);
ir_partial=ir(1:Ntx)/max(abs(ir));

%%--- BW : for figures
switch wav
    case 1, BW=f_IF/1000+1/PW;   % [MHz] f_IF + (null-to-null BW of pulse)
    case 2, BW=1/PW;        % [MHz] null-to-null BW of rectangular pulse
    case 3, BW=2*f_IF/1000;       % [MHz] spreading BW of linear chirp
end

%%--- flags
rflag=0; % if outside of maxmum range, flag=1
signal_miss=1; % if target missing, flag==1
only_ang_miss=0; % if angle missing, flag==1
jflag=zeros(1,num); % jamming flag

%%--- initialization
ang_miss_cnt=0; % for re-tracking after angle miss
t_sum=0; % for velocity jamming
R_sqerr_sum=0; V_sqerr_sum=0; % rms errors
Az_sqerr_sum=0; Ev_sqerr_sum=0; % rms errors

%%--- target moving position 
target=Fn_TargetMove(ip_target,vel,ad_target,num,t_intv);


%% ====== Start Simulations =================
for k=1:num

%% *** Received Signal ***

%%--- Generate (k)th Target Position
[R, azim, elev]=Fn_RangeAzEl(target(k,:),h_target-h_radar);
% R [km] : target range
% azim [deg] : target azimuth angle
% elev [deg] : target elevation angle
azim_rad=azim*pi/180; % [rad] target azimuth angle
elev_rad=elev*pi/180; % [rad] target elevation angle
if R>max_range, rflag=1; else, rflag=0; end

%%--- Doppler Frequency
[fd, direction]=Fn_DopplerFreq(target(k,:),h_radar,h_target,ad_target,vel,f_RF);
% fd [Hz] : Doppler frequency 
% direction : direction vector
%%--- Pulse Width (spreading)
psi_g=atan(2*h_radar/(R*cos(elev_rad))+tan(elev_rad)); % [rad] grazing angle
G=Fn_Reflection(f_RF, psi_g, ground_type, rainfall, temper-(6.7*h_target)); % G=[Gv,Gh]
rho3=Fn_PulseSpread(elev_rad,psi_g,R,G,theta3); % [micro sec] delay spread
rho3=round(rho3/2/tt)*2*tt; % [micro sec] even muliple of tt

%-- Pulse Power (two way)
PA_rcv=Fn_PulseAmp(PP-30+Gt+Gr,PRI,PW,R,f_RF,RCS,h_target,rainfall);  
% % %--- for Combined Simulator (5F Lab)
% % f_shift=f_IF+fd*10^(-6);
% % PW_rcv=PW+rho_half*2; % two-way
% % PP_rcv=P+Gt+Gr-Lfs; % received power, no consideration of weather, dB-W

%%--- Waveform
txwave=t; Ntx=Nt;
tr_wav=Fn_Waveform(wav,f_MF,0,PW,txwave,rho3,0);
pulse_spread=0;
if pulse_spread==0
    rxwave=txwave; 
    Nrx=Ntx;
    rcv_wav=Fn_Waveform(wav,f_MF,fd*10^(-6),PW,txwave,rho3,1);
    rcv_wav_IF=Fn_Waveform(wav,f_IF/1000,fd*10^(-6),PW,txwave,rho3,1);    
else
    rxwave=0:tt:PW+rho3-tt; % [micro sec] pulse (waveform) time
    Nrx=length(rxwave);
	rcv_wav=Fn_Waveform(wav,f_MF,fd*10^(-6),PW,rxwave,rho3,2);
	rcv_wav_IF=Fn_Waveform(wav,f_IF/1000,fd*10^(-6),PW,rxwave,rho3,2); 
end

%%--- Received signal (PRI 길이에서)
t_target=R/c*10*2; % [micro sec] target time (multiple of tt)
if t_target<PRI
    start_rcv=0; % [micro sec]
    end_rcv=PRI; % [micro sec]
else
    t_target=t_target-PRI;
    start_rcv=PRI; % [micro sec]
    end_rcv=PRI*2; % [micro sec]
end
Ntt=round(t_target/tt);
rcv_pulse=PA_rcv*rcv_wav;
rcv_IFpulse=PA_rcv*rcv_wav_IF;
Nrcvp=length(rcv_pulse);
t_rcv=(0:Nrcvp-1)*tt+start_rcv;

%%%% Doppler Freq View %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tx_spec=abs(fft(rcv_IFpulse,FBank));
% tx_spec=tx_spec/max(tx_spec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% *** Jamming ***
if type_jam>0
    dur_jam=dur_jamming/t_intv;
    delay_jam=delay_jamming/t_intv;
    if R<=max_range && jflag(k)==0 % jflag setting, no jamming
        total_jam=delay_jam+dur_jam;
        end_jam=min(k+total_jam-1,num);
        jflag(k:end_jam)=total_jam:-1:(total_jam-end_jam+k);
        t_sum=0;
    elseif jflag(k)<=dur_jam && jflag(k)>0 % do jamming
        Pj=10^((PP-30+Gt)/10)*RCS/((BT_range*10^3)^2*(4*pi)*10^(Gj/10)); % [Watt]
        switch type_jam
            case 11 % barrage noise
                amp_jam=Pj/(wc2-wc1);
                barr_ir=amp_jam*sinc(2*(wc2-wc1)*(-150:150)*tt);
                tmp=conv(randn(1,Nrcvp),barr_ir);
                rcv_pulse=rcv_pulse+tmp(151:Nrcvp+150);
                rcv_IFpulse=rcv_IFpulse+tmp(151:Nrcvp+150);
            case 12 % spot noise
                amp_jam=sqrt(Pj/spot_num*2)*(PW/PRI)*tt; % [Watt] spot noise power
                spot_freq=f_IF/1000+((1:spot_num)-round(spot_num/2))*0.03; % [MHz
                for ns=1:spot_num
                    rcv_pulse=rcv_pulse+cos(2*pi*spot_freq(ns)*t_rcv)*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+cos(2*pi*spot_freq(ns)*t_rcv)*amp_jam;
                end
            case 21 % RGPO
                ;
            case 22
                ;
            case 31
                v_jam=vel-acc_jam*(dur_jam-jflag(k))*t_intv; % [m/s] 
                f_jam=(f_RF+2*direction*v_jam/lambda*10^-6)/(1+vel*direction/(c*10^8)); % [MHz]
                fd_jam=f_jam*vel*direction/c*10^(-2); % [Hz] Doppler freq of f_jam
                fd_radar=f_jam+fd_jam*10^(-6)-f_RF; % [MHz] Dopper freq estimated by radar
%                 [vel, v_jam, c*(10^2)*fd_radar/(2*direction*f_RF)]
                jrcv_wav=Fn_Waveform(wav,f_MF,fd_radar,PW,txwave,rho3,1);
                jrcv_wav_IF=Fn_Waveform(wav,f_IF/1000,fd_radar,PW,txwave,rho3,1);
                t_jam=t_target;
                amp_jam=sqrt(Pj*2)*(PW/PRI)*tt;
                rcv_pulse=rcv_pulse+jrcv_wav*amp_jam;
                rcv_IFpulse=rcv_IFpulse+jrcv_wav_IF*amp_jam;                
            case 32 %% VGPI
                v_jam=vel+acc_jam*(dur_jam-jflag(k))*t_intv; % [m/s] 
                f_jam=(f_RF+2*direction*v_jam/lambda*10^-6)/(1+vel*direction/(c*10^8)); % [MHz]
                fd_jam=f_jam*vel*direction/c*10^(-2); % [Hz] Doppler freq of f_jam
                fd_radar=f_jam+fd_jam*10^(-6)-f_RF; % [MHz] Dopper freq estimated by radar
%                 [vel, v_jam, c*(10^2)*fd_radar/(2*direction*f_RF)]
                jrcv_wav=Fn_Waveform(wav,f_MF,fd_radar,PW,txwave,rho3,1);
                jrcv_wav_IF=Fn_Waveform(wav,f_IF/1000,fd_radar,PW,txwave,rho3,1);
                t_jam=t_target;
                amp_jam=sqrt(Pj*2)*(PW/PRI)*tt;
                rcv_pulse=rcv_pulse+jrcv_wav*amp_jam;
                rcv_IFpulse=rcv_IFpulse+jrcv_wav_IF*amp_jam;
            case 61 %% VGPO 
                v_jam=vel-acc_jam*(dur_jam-jflag(k))*t_intv; % [m/s] 
                f_jam=(f_RF+2*direction*v_jam/lambda*10^-6)/(1+vel*direction/(c*10^8)); % [MHz]
                fd_jam=f_jam*vel*direction/c*10^(-2); % [Hz] Doppler freq of f_jam
                fd_radar=f_jam+fd_jam*10^-6-f_RF; % [MHz] Dopper freq estimated by radar
                jrcv_wav=Fn_Waveform(wav,f_MF,fd_radar,PW,txwave,rho3,1);
                jrcv_wav_IF=Fn_Waveform(wav,f_IF/1000,fd_radar,PW,txwave,rho3,1);
                t_sum=t_sum+(dur_jam-jflag(k))*acc_jam*2*direction/c*10^-2;
                t_jam=t_target-t_sum;
                % t_jam=t_target+(dur_jam-jflag(k))^2*acceleration*direction/c*10^-2; % integration
                amp_jam=sqrt(Pj*2)*(PW/PRI)*tt;                
                rcv_pulse=rcv_pulse+jrcv_wav*amp_jam;
                rcv_IFpulse=rcv_IFpulse+jrcv_wav_IF*amp_jam;
            case 62 %% VGPI
                v_jam=vel+acc_jam*(dur_jam-jflag(k))*t_intv; % [m/s] 
                f_jam=(f_RF+2*direction*v_jam/lambda*10^-6)/(1+vel*direction/(c*10^8)); % [MHz]
                fd_jam=f_jam*vel*direction/c*10^(-2); % [Hz] Doppler freq of f_jam
                fd_radar=f_jam+fd_jam*10^(-6)-f_RF; % [MHz] Dopper freq estimated by radar
                jrcv_wav=Fn_Waveform(wav,f_MF,fd_radar,PW,txwave,rho3,1);
                jrcv_wav_IF=Fn_Waveform(wav,f_IF/1000,fd_radar,PW,txwave,rho3,1);
                t_sum=t_sum+(dur_jam-jflag(k))*acc_jam*2*direction/c*10^-2;
                t_jam=t_target+t_sum;
                % t_jam=t_target+(dur_jam-jflag(k))^2*acceleration*direction/c*10^-2; % integration
                amp_jam=sqrt(Pj*2)*(PW/PRI)*tt;
                rcv_pulse=rcv_pulse+jrcv_wav*amp_jam;
                rcv_IFpulse=rcv_IFpulse+jrcv_wav_IF*amp_jam;
        end
    end
end

%% *** Radar Receiver ***
% %%--- Anttena Gain
% % % G=g0*exp(-log(2)*(2*the/theta3).^2) % Gaussian beam
% if k==1
%     az_beam=azim-rand(1,1)*2; % initial beam center
%     el_beam=elev-rand(1,1)*2; % initial beam center
% elseif ang_miss_cnt>5 % start re-tracking
%     az_beam=azim-rand(1,1)*2; % initial beam center
%     el_beam=elev-rand(1,1)*2; % initial beam center    
% end
% if azim>315 && az_beam<45 %%% 0 <--> 360 
%     az_beam=az_beam+360; 
% elseif azim<45 && az_beam>315
%     az_beam=az_beam-360;
% end
% azi_angle=azim-az_beam; % target azi angle from beam center, deg
% elv_angle=elev-el_beam; % target elv angle from beam center, deg
% km=2*sqrt(2)*log(2)*theta_s/theta3; % antenna constant
% V=Fn_AnttenaGain(g0, theta3, theta_s, azi_angle, elv_angle);
% % if type_radar==2 || type_radar==3
%     rcv_signal=( V(1)*(rcv_pulse+randn(1,Nrcvp)*sigma) ...
%         + V(2)*(rcv_pulse+randn(1,Nrcvp)*sigma) ...
%         + V(3)*(rcv_pulse+randn(1,Nrcvp)*sigma) ...
%         + V(4)*(rcv_pulse+randn(1,Nrcvp)*sigma) )/sum(V);    
%     rcv_IFsignal=( V(1)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) ...
%         + V(2)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) ...
%         + V(3)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) ...
%         + V(4)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) )/sum(V);  
% % else
    rcv_signal=g0*rcv_pulse+randn(1,Nrcvp)*sigma;
    rcv_IFsignal=g0*rcv_IFpulse+randn(1,Nrcvp)*sigma;
% end

%%--- Matched Filter (Range Estimation)
Csr=[];
Ntx=PRI/4/tt;
for tau=1:Nrcvp-Ntx+1
    rcv_window = rcv_IFsignal(tau:tau+Ntx-1);
    Csr(tau)=Fn_CorrAbs(rcv_wav_IF(1:Ntx),rcv_window);
end
[mmcsr,idx]=max(Csr);
est_R=(idx*tt+start_rcv)*c/10*0.5; % [km]
if mmcsr>Vt 
    if signal_miss==1 
        detect_R=est_R;
    end % first detecting position 
    signal_miss=0;
else
    signal_miss=1;
    rcv_signal=randn(1,Nrcvp)*sigma;
    rcv_IFsignal=rcv_signal;
end

% %%--- Monopulse Diff and Sum (Angle Estimation)
% if signal_miss==0
%     delta_az=(V(1)-V(2))-(V(3)-V(4));
%     delta_el=(V(1)+V(2))-(V(3)+V(4)); 
%     hat_azi=delta_az*theta3/(4*km); % estimation of azi_angle
%     hat_elv=delta_el*theta3/(4*km); % estimation of elv_angle
%     if type_radar==2 || type_radar==3
%         if abs((az_beam+hat_azi)-azim) <= theta3 ...
%                 || (az_beam>350 && abs((az_beam+hat_azi-360)-azim) <= theta3)
%             only_ang_miss=0;
%         else
%             only_ang_miss=1; % angle tracking missing
%             rcv_signal=randn(1,Nrcvp)*sigma;
%             rcv_IFsignal=rcv_signal;
%             ang_miss_cnt=ang_miss_cnt+1; % for re-tracking after angle miss
%         end
%     end
% end

%%--- IF and Limiter / AGC
switch gain_control
    case 1 %% Limiter
        gc_signal=rcv_IFsignal/max(rcv_IFsignal);
    case 2
        %%-- IF
        real_op_range=Vt*10^(op_range/10); % operating range
        if max(rcv_signal)<real_op_range
            IF_signal=10^(G_IF/10)*real(conv(ir_partial,rcv_IFsignal));
        else
            IF_signal=rcv_IFsignal;
        end
%         IF_signal=real(IF_signal(1:Nrcvp)*tt);
        %%-- AGC
        switch agc_dyn
            case 1, dyn_range=[-50 50]; % [dB]
            case 2, dyn_range=[-30 30]; % [dB]
        end
        peak_dB=10*log10(mmcsr)+G_IF;
        if peak_dB>dyn_range(1) && peak_dB<dyn_range(2)
            if k<(agc_delay/t_intv)
                gc_signal=IF_signal/max(IF_signal)*k/(agc_delay/t_intv);
            else
                gc_signal=IF_signal/max(IF_signal);
            end
        else % 범위 안에 들지 않으면 크기 맞춰지지 않음.
            gc_signal=IF_signal;
        end
end

% %%-- Range Tracking
% if (type_radar==2 || type_radar==3) && signal_miss==0
%     if k==1, r_gate=idx:idx+Ntx-1; end
%     e_gate=r_gate(1):r_gate(Ntx/2);
%     l_gate=r_gate(Ntx/2+1):r_gate(end);
%     eg_sum=sum(abs(gc_signal(e_gate)));
%     lg_sum=sum(abs(gc_signal(l_gate)));
%     r_gate=r_gate+round(lg_sum-eg_sum);
%     est_Rg=(r_gate(1)*tt+start_rcv)*c/10*0.5; % [km]
%     if abs(est_R-est_Rg)>0.5, est_R=est_Rg; end
% end

%%--- BPF Bank (Frequency Estimation)
rcv_spec=abs(fft(gc_signal,FBank)); % full signal
[mmspec,idx_f]=max(rcv_spec);
est_fd=10^6*((idx_f-1)/(tt*FBank))-f_IF*1000; % [Hz] estimated Doppler freq
est_vel=c*(10^2)*est_fd/(2*direction*f_RF); % [m/s] estimated velocity


switch ccc
case 0

    %%--- figures for rcv_signal, Csr, Vt, IF_signal, agc_signal, etc 
clf
subplot(321)   
    plot(t_rcv, rcv_signal); title('rcv signal');
    rmm=max(abs(rcv_signal)); axis([start_rcv end_rcv -rmm rmm]); 
subplot(322) 
    if start_rcv==0 
        plot((1:length(Csr))*tt-tt, Csr), hold on,         
    else
        plot(PRI+(1:length(Csr))*tt-tt, Csr), hold on, 
    end
    plot([start_rcv end_rcv],[Vt Vt],'r'), title('correlation');
    rmm=Vt*10; axis([start_rcv end_rcv 0 rmm]);    
%     text(start_rcv+10, 4*10^-3, ['R = ',num2str(R),' km']);
%     text(start_rcv+10, 3.5*10^-3, ['R_{est} = ',num2str(est_R),' km']);
subplot(323)
    plot(t_rcv, rcv_IFsignal); title('rcv IFsignal');
    rmm=max(abs(rcv_IFsignal)); axis([start_rcv end_rcv -rmm rmm]); 
subplot(324)
    plot(t_rcv, real(gc_signal(1:Nrcvp))); title('gain controlled signal')
    rmm=max(abs(gc_signal)); axis([start_rcv end_rcv -rmm rmm])      
subplot(325)
    rcv_spec=abs(fft(rcv_signal,FBank));
    rcv_spec=rcv_spec/max(rcv_spec);
    semilogy((0:length(rcv_spec)-1)/(tt*FBank),rcv_spec); grid on;   
    axis([0 min(f_MF*5,Fs/2) 10^(-3) 2]);
    xlabel('frequency (GHz)'), title('rcv spectrum');
subplot(326)  
    rcv_spec=abs(fft(gc_signal,FBank));
    rcv_spec=rcv_spec/max(rcv_spec);
    semilogy((0:length(rcv_spec)-1)/(tt*FBank),rcv_spec); grid on;   
    axis([0 min(f_IF*5/1000,Fs/2) 10^(-4) 2]);
    xlabel('frequency (MHz)'), title('gain controlled spectrum');
%     text(f_IF*2/1000, 10^(0), ['v = ',num2str(vel),' m/s']);
%     text(f_IF*2/1000, 10^(-1)*4, ['v_{est} = ',num2str(est_vel),' m/s']);
%     rcv_spec=abs(fft(rcv_wav_IF,FBank));
%     rcv_spec=rcv_spec/max(rcv_spec);
%     hold on; semilogy((0:length(rcv_spec)-1)/(tt*FBank),rcv_spec); 
pause(t_intv/proc_speed);

case 1
    
%% *** Plotting ***
est=Fn_TargetPoSearch(target(k,:),h_target,theta3,est_R);
x_est=est(1); y_est=est(2);
[r2, a2, e2]=Fn_RangeAzEl([x_est,y_est],h_target);
est_bin_az=floor(a2/theta3)*theta3;
est_bin_el=floor(e2/theta3)*theta3;
n_est_bin=floor(est_bin_az/theta3)+1; % n-th bin

figure(1), clf;
if mod(k,Tscan)==n_est_bin
    
subplot(221) %%%%%%%%%%%%
    Fn_PlotSearchScope(max_range,target(k,:),BT_range);
    if jflag(k)>0 && jflag(k)<=dur_jam && type_jam>0
        plot(target(k,1),target(k,2),'go','MarkerFaceColor','g',...
            'MarkerSize',15);
        plot([target(k,1) -max_range],[target(k,2) -max_range],':','Color','green')
        text(-max_range-10, -max_range-5,'jamming','Color','green');
        plot(target(k,1),target(k,2),'rx','LineWidth',2);
    end
    if rflag==0 && signal_miss==0 % && only_ang_miss==0
        [r2, a2, e2]=Fn_RangeAzEl([x_est,y_est],h_target);
        est_bin=floor(a2/theta3)*theta3;
        est_bin_el=floor(e2/theta3)*theta3;
        plot([0 max_range*cos(est_bin*pi/180)],...
            [0 max_range*sin(est_bin*pi/180)],'Color','blue');
        plot([0 max_range*cos((est_bin+theta3)*pi/180)],...
            [0 max_range*sin((est_bin+theta3)*pi/180)],'Color','blue');
        if est_bin<0
            est_bin=est_bin+360;
        elseif est_bin>360
            est_bin=est_bin-360;
        end
        if est_bin<180
            if est_bin+theta3<=180 %%% mod(180,theta3)~=0 일때 필요함
                cx=max_range*cos((est_bin+theta3)*pi/180):0.01:max_range*cos(est_bin*pi/180);
                plot(cx,sqrt(max_range^2-cx.^2),'Color','blue','LineWidth',3);
            else
                cx=-max_range:0.01:max_range*cos(est_bin*pi/180);
                plot(cx,sqrt(max_range^2-cx.^2),'Color','blue','LineWidth',3);
                cx=-max_range:0.01:max_range*cos((est_bin+theta3)*pi/180);
                plot(cx,-sqrt(max_range^2-cx.^2),'Color','blue','LineWidth',3);
            end
        else % est_bin>=180
            if est_bin+theta3<=360 %%% mod(180,theta3)~=0 일때 필요함
                cx=max_range*cos(est_bin*pi/180):0.01:max_range*cos((est_bin+theta3)*pi/180);
                plot(cx,-sqrt(max_range^2-cx.^2),'Color','blue','LineWidth',3);
            else
                cx=max_range*cos(est_bin*pi/180):0.01:max_range;
                plot(cx,-sqrt(max_range^2-cx.^2),'Color','blue','LineWidth',3);
                cx=max_range*cos((est_bin+theta3)*pi/180):0.01:max_range;
                plot(cx,sqrt(max_range^2-cx.^2),'Color','blue','LineWidth',3);
            end
        end
    else
        plot([0 max_range*cos((k-1)*theta3*pi/180)],...
            [0 max_range*sin((k-1)*theta3*pi/180)],...
            'Color','black','LineWidth',2);
        plot([0 max_range*cos(k*theta3*pi/180)],...
            [0 max_range*sin(k*theta3*pi/180)],...
            'Color','black','LineWidth',2);
    end
%     text(-max_range/2,-max_range-10,'CW Doppler',...
%         'Color','black','FontSize',14);

subplot(222) %%%%%%%%%%%%
    if rflag==0 && signal_miss==0 % && only_ang_miss==0
        Fn_PlotAngScope(theta3,est_bin,est_bin_el,signal_miss);
        plot(azim,elev,'rx','LineWidth',2),
%         text(est_bin+theta3/4,est_bin_el-theta3/4, ...
%             'CW Dopper','Color','black','FontSize',14);
    else
        Fn_PlotAngScope(theta3,0,0,1);
%         text(-theta3/4,-theta3/5*4, ...
%             'CW Dopper','Color','black','FontSize',14);
    end

graph_signal=gc_signal(1:Nrcvp);
spec_signal=gc_signal;
rmm=1.2;
% rmm=max(abs(graph_signal))*1.2; 
subplot(223),
    if rflag==1 || only_ang_miss==1
        plot(t_rcv, graph_signal), grid on; hold on;
        axis([start_rcv end_rcv -rmm, rmm]);
    else
        plot(t_rcv, graph_signal), grid on; hold on;
        axis([start_rcv end_rcv -rmm rmm]),
    end
    title('Received Signal'), xlabel('time (\mus)');
    
subplot(224),    
%%%% Doppler frequency View %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     semilogy((0:length(tx_spec)-1)/(tt*FBank),tx_spec); hold on; 
%     axis([0.495 0.515 10^(-2) 1.5]);  
%         text(0.508, 10^(0), ['v = ',num2str(vel),' m/s']);
%         text(0.508, 10^(-1)*7, ['v_{est} = ',num2str(est_vel),' m/s']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rcv_spec=abs(fft(spec_signal,FBank));
    rcv_spec=rcv_spec/max(rcv_spec);
    semilogy((0:length(rcv_spec)-1)/(tt*FBank),rcv_spec); grid on;   
    axis([0 min(5,Fs/2) 10^(-4) 2]);
    xlabel('frequency (MHz)'), title('Normalized Spectrum');
    if rflag==0 && signal_miss==0 && only_ang_miss==0 
        text(BW*2, 10^(0), ['v = ',num2str(vel),' m/s']);
        text(BW*2, 10^(-1)*5, ['v_{est} = ',num2str(est_vel),' m/s']);
    end
    

% %%%%%%%%%% for presentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2), clf
%         Fn_PlotSearchScope(max_range,target(k,:),BT_range);
%         if jflag(k)>0 && jflag(k)<=dur_jam && type_jam>0 %%% jamming
%             plot(target(k,1),target(k,2),'go','MarkerFaceColor','g',...
%                 'MarkerSize',15);  
%             plot([target(k,1) -max_range],[target(k,2) -max_range],':','Color','green')
%             text(-max_range-10, -max_range-5,'jamming','Color','green');
%             plot(target(k,1),target(k,2),'rx','LineWidth',2);
%         end
%         if rflag==0 && signal_miss==0 && only_ang_miss==0
%             [r2, a2, e2]=Fn_RangeAzEl([x_est,y_est],h_target);
%             est_bin=floor(a2/theta3)*theta3;
%             plot([0 max_range*cos(est_bin*pi/180)],...
%                 [0 max_range*sin(est_bin*pi/180)],'Color','blue');
%             plot([0 max_range*cos((est_bin+theta3)*pi/180)],...
%                 [0 max_range*sin((est_bin+theta3)*pi/180)],'Color','blue');    
%             if est_bin<180
%                 if est_bin+theta3<=180 %%% mod(180,theta3)~=0 일때 필요함
%                     cx=r2*cos((est_bin+theta3)*pi/180):0.01:r2*cos(est_bin*pi/180);
%                     plot(cx,sqrt(r2^2-cx.^2),'Color','blue','LineWidth',3);
%                 else
%                     cx=-r2:0.01:r2*cos(est_bin*pi/180);
%                     plot(cx,sqrt(r2^2-cx.^2),'Color','blue','LineWidth',3);
%                     cx=-r2:0.01:r2*cos((est_bin+theta3)*pi/180);
%                     plot(cx,-sqrt(r2^2-cx.^2),'Color','blue','LineWidth',3);
%                 end
%             else % est_bin>=180
%                 if est_bin+theta3<=360 %%% mod(180,theta3)~=0 일때 필요함
%                     cx=r2*cos(est_bin*pi/180):0.01:r2*cos((est_bin+theta3)*pi/180);
%                     plot(cx,-sqrt(r2^2-cx.^2),'Color','blue','LineWidth',3);
%                 else
%                     cx=r2*cos(est_bin*pi/180):0.01:r2;
%                     plot(cx,-sqrt(r2^2-cx.^2),'Color','blue','LineWidth',3);
%                     cx=r2*cos((est_bin+theta3)*pi/180):0.01:r2;
%                     plot(cx,sqrt(r2^2-cx.^2),'Color','blue','LineWidth',3);
%                 end
%             end
%         end
%         text(-max_range-10,-max_range-10,'Search',...
%             'Color','black','FontSize',14);
%         plot(target(k,1),target(k,2),'rx','LineWidth',2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% *** Search/Tracking Results ***
if rflag==1
    'outside of maximum range'
elseif signal_miss==1
    'target missing'
elseif only_ang_miss==1
    'angle missing'
else
    [vel,est_vel]
    V_sqerr_sum=V_sqerr_sum+(vel-est_vel)^2;
end
        
else %% bin
rcv_IFsignal=randn(1,Nrcvp)*sigma;
gc_signal=rcv_IFsignal/max(rcv_IFsignal);

subplot(221) %%%%%%%%%%%%
        Fn_PlotSearchScope(max_range,target(k,:),BT_range);
        if jflag(k)>0 && jflag(k)<=dur_jam && type_jam>0 %%% jamming
            plot(target(k,1),target(k,2),'go','MarkerFaceColor','g',...
                'MarkerSize',15);  
            plot([target(k,1) -max_range],[target(k,2) -max_range],':','Color','green')
            text(-max_range-10, -max_range-5,'jamming','Color','green');
            plot(target(k,1),target(k,2),'rx','LineWidth',2);
        end
        plot([0 max_range*cos((k-1)*theta3*pi/180)],...
            [0 max_range*sin((k-1)*theta3*pi/180)],...
            'Color','black','LineWidth',2);
        plot([0 max_range*cos(k*theta3*pi/180)],...
            [0 max_range*sin(k*theta3*pi/180)],...
            'Color','black','LineWidth',2);
%         text(-max_range/2,-max_range-10,'Search','Color','black','FontSize',14);
%         plot(target(k,1),target(k,2),'rx','LineWidth',2);

subplot(222) %%%%%%%%%%%%
	Fn_PlotAngScope(theta3,0,0,1);

graph_signal=gc_signal(1:Nrcvp);
spec_signal=gc_signal;
rmm=1.2;
% rmm=max(abs(graph_signal))*1.2; 
subplot(223),
    if rflag==1 || only_ang_miss==1
        plot(t_rcv, graph_signal), grid on; hold on;
        axis([start_rcv end_rcv -rmm, rmm]);
    else
        plot(t_rcv, graph_signal), grid on; hold on;
        axis([start_rcv end_rcv -rmm rmm]),
    end
    title('Received Signal'), xlabel('time (\mus)');
    if rflag==0 && signal_miss==0 && only_ang_miss==0 ...
            && (type_radar==2 || type_radar==3)
        text(start_rcv+10, rmm*0.8, ['R = ',num2str(R),' km']);
        text(start_rcv+10, rmm*0.6, ['R_{est} = ',num2str(est_R),' km']);
%         if (r_gate(1)+Ntx-1)>Nrcvp
%             plot(t_rcv, [zeros(1,r_gate(1)-1),ones(1,Nrcvp-r_gate(1)+1)]);
%         else
%             plot(t_rcv, [zeros(1,r_gate(1)-1),ones(1,Ntx), ...
%                 zeros(1,Nrcvp-r_gate(1)-Ntx+1)]);
%         end
    end
    
subplot(224),    
%%%% Doppler frequency View %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     semilogy((0:length(tx_spec)-1)/(tt*FBank),tx_spec); hold on; 
%     axis([0.495 0.515 10^(-2) 1.5]);  
%         text(0.508, 10^(0), ['v = ',num2str(vel),' m/s']);
%         text(0.508, 10^(-1)*7, ['v_{est} = ',num2str(est_vel),' m/s']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rcv_spec=abs(fft(spec_signal,FBank));
    rcv_spec=rcv_spec/max(rcv_spec);
    semilogy((0:length(rcv_spec)-1)/(tt*FBank),rcv_spec); grid on;   
    axis([0 min(5,Fs/2) 10^(-4) 2]);
    xlabel('frequency (MHz)'), title('Normalized Spectrum');
    
end %% bin

end %% ccc

pause(t_intv/proc_speed);
end %% end of num
V_rms=sqrt(V_sqerr_sum/num)     % rms error of estimated vel [m/s]