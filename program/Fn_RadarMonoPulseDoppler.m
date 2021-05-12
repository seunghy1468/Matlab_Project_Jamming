function Fn_RadarMonoPulseDoppler()

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
if Ntt+Nrx>Nt
    rcv_pulse=PA_rcv*[zeros(1,Ntt-1),rcv_wav];
    Nrcvp=length(rcv_pulse);
    %         rcv_signal=rcv_pulse+randn(1,Nrcv)*sigma;
    rcv_IFpulse=PA_rcv*[zeros(1,Ntt-1),rcv_wav_IF];
else
    rcv_pulse=PA_rcv*[zeros(1,Ntt-1),rcv_wav,zeros(1,Nt-(Ntt+Nrx-1))];
    Nrcvp=length(rcv_pulse);
    rcv_IFpulse=PA_rcv*[zeros(1,Ntt-1),rcv_wav_IF,zeros(1,Nt-(Ntt+Nrx-1))];
end
t_rcv=(0:Nrcvp-1)*tt+start_rcv;

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
                t_jam=t_target+(dur_jam-jflag(k))*t_intv*vel_jam/c*10^(-2)*2;
                amp_jam=sqrt(Pj*2)*(PW/PRI)*tt;
                if (Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)<1
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                else
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                end
                
            case 22
                t_jam=t_target-(dur_jam-jflag(k))*t_intv*vel_jam/c*10^(-2)*2;
                amp_jam=sqrt(Pj*2)*(PW/PRI)*tt;
                if (Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)<1
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                else
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                end
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
                if (Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)<1
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                else
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                end
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
                if (Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)<1
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                else
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                end
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
                if (Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)<1
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                else
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                end
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
                if (Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)<1
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF(1:Nrcvp-(round((t_jam-start_rcv)/tt)-1))]*amp_jam;
                else
                    rcv_pulse=rcv_pulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                    rcv_IFpulse=rcv_IFpulse+[zeros(1,round((t_jam-start_rcv)/tt)-1),...
                        rcv_wav_IF,zeros(1,Nrcvp-(round((t_jam-start_rcv)/tt)-1)-Nrx)]*amp_jam;
                end
        end
    end
end

%% *** Radar Receiver ***
%%--- Anttena Gain
% % G=g0*exp(-log(2)*(2*the/theta3).^2) % Gaussian beam
if k==1
    az_beam=azim-rand(1,1)*2; % initial beam center
    el_beam=elev-rand(1,1)*2; % initial beam center
elseif ang_miss_cnt>5 % start re-tracking
    az_beam=azim-rand(1,1)*2; % initial beam center
    el_beam=elev-rand(1,1)*2; % initial beam center    
end
if azim>315 && az_beam<45 %%% 0 <--> 360 
    az_beam=az_beam+360; 
elseif azim<45 && az_beam>315
    az_beam=az_beam-360;
end
azi_angle=azim-az_beam; % target azi angle from beam center, deg
elv_angle=elev-el_beam; % target elv angle from beam center, deg
km=2*sqrt(2)*log(2)*theta_s/theta3; % antenna constant
V=Fn_AnttenaGain(g0, theta3, theta_s, azi_angle, elv_angle);
rcv_signal=( V(1)*(rcv_pulse+randn(1,Nrcvp)*sigma) ...
    + V(2)*(rcv_pulse+randn(1,Nrcvp)*sigma) ...
    + V(3)*(rcv_pulse+randn(1,Nrcvp)*sigma) ...
    + V(4)*(rcv_pulse+randn(1,Nrcvp)*sigma) )/sum(V);
rcv_IFsignal=( V(1)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) ...
    + V(2)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) ...
    + V(3)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) ...
    + V(4)*(rcv_IFpulse+randn(1,Nrcvp)*sigma) )/sum(V);

%%--- Matched Filter (Range Estimation)
Csr=[];
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

%%--- Monopulse Diff and Sum (Angle Estimation)
if signal_miss==0
    delta_az=(V(1)-V(2))-(V(3)-V(4));
    delta_el=(V(1)+V(2))-(V(3)+V(4));
    hat_azi=delta_az*theta3/(4*km); % estimation of azi_angle
    hat_elv=delta_el*theta3/(4*km); % estimation of elv_angle
    if abs((az_beam+hat_azi)-azim) <= theta3 ...
            || (az_beam>350 && abs((az_beam+hat_azi-360)-azim) <= theta3)
        only_ang_miss=0;
    else
        only_ang_miss=1; % angle tracking missing
        rcv_signal=randn(1,Nrcvp)*sigma;
        rcv_IFsignal=rcv_signal;
        ang_miss_cnt=ang_miss_cnt+1; % for re-tracking after angle miss
    end
end

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

%%-- Range Tracking
if (type_radar==2 || type_radar==3) && signal_miss==0
    if k==1, r_gate=idx:idx+Ntx-1; end
    e_gate=r_gate(1):r_gate(Ntx/2);
    l_gate=r_gate(Ntx/2+1):r_gate(end);
    eg_sum=sum(abs(gc_signal(e_gate)));
    lg_sum=sum(abs(gc_signal(l_gate)));
    r_gate=r_gate+round(lg_sum-eg_sum);
    est_Rg=(r_gate(1)*tt+start_rcv)*c/10*0.5; % [km]
    if abs(est_R-est_Rg)>0.5, est_R=est_Rg; end
end

%%--- BPF Bank (Frequency Estimation)
rcv_spec=abs(fft(gc_signal(Ntt:Ntt+Ntx-1),FBank)); % pulse signal
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
figure(1), clf;

switch type_jam
    case 11, str_jam='Barrage Noise Jamming';
    case 12, str_jam='Spot Noise Jamming';
    case 21, str_jam='Range Gate Pull Off'; % 멀어짐
    case 22, str_jam='Range Gate Pull In';  % 가까워짐
    case 31, str_jam='Velocity Gate Pull Off'; % 빨라짐
    case 32, str_jam='Velocity Gate Pull In';  % 느려짐
    case 62, str_jam='RGPO + Velocity Deception'; % 빠르게 멀어짐
    case 61, str_jam='RGPI + Velocity Deception'; % 느리게 멀어짐
    otherwise, str_jam='No Jamming';
end
% Position: [a b c d]
% 왼쪽 아래 꼭지점 위치: 
% - 그림의 왼쪽에서부터 가로 길이의 a 비율 만큼 떨어진 곳
% - 그림의 아래쪽에서부터 세로 길이의 b 비율 만큼 떨어진 곳
% 오른쪽 위 꼭지점 위치:
% - 왼쪽 아래 꼭지점(a)에서부터 가로 길이의 c 비율 만큼 떨어진 곳
% - 왼쪽 아래 꼭지점(b)에서부터 세로 길이의 d 비율 만큼 떨어진 곳
p1 = uipanel('Position',[.66 .82 .25 .15]);
uicontrol(p1,'Style','text','Units','normalized',...
               'String','<시뮬레이션 환경>','FontSize',12,...
               'Position',[.15 .6 .7 .3]);
uicontrol(p1,'Style','text','Units','normalized',...
               'String','위협체: Monopulse Doppler Radar','FontSize',10,...
               'Position',[.1 .35 .85 .2],'HorizontalAlignment','left');
uicontrol(p1,'Style','text','Units','normalized',...
               'String',['재밍: ',str_jam],'FontSize',10,...
               'Position',[.1 .1 .8 .2],'HorizontalAlignment','left');  
p2 = uipanel('Position',[.66 .25 .25 .55]);
uicontrol(p2,'Style','text','Units','normalized',...
               'String','<추적 상황>','FontSize',12,...
               'Position',[.15 .87 .7 .1]);
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['거리: ',num2str(R),' km'],'FontSize',10,...
               'Position',[.1 .82 .8 .05],'HorizontalAlignment','left');   
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['방위각: ',num2str(azim),' deg'],'FontSize',10,...
               'Position',[.1 .60 .8 .05],'HorizontalAlignment','left');  
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['고도각: ',num2str(elev),' deg'],'FontSize',10,...
               'Position',[.1 .38 .8 .05],'HorizontalAlignment','left');   
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['속도: ',num2str(vel),' m/s'],'FontSize',10,...
               'Position',[.1 .16 .8 .05],'HorizontalAlignment','left');

subplot(231) %%%%%%%%%%%%
    Fn_PlotTrackScope(max_range,target(k,:),BT_range);
    if jflag(k)>0 && jflag(k)<=dur_jam && type_jam>0
        plot(target(k,1),target(k,2),'go','MarkerFaceColor','g',...
            'MarkerSize',15);
        plot([target(k,1) -max_range],[target(k,2) -max_range],':','Color','green')
        text(-max_range-10, -max_range-5,'jamming','Color','green');
        plot(target(k,1),target(k,2),'rx','LineWidth',2);
    end
    if rflag==0 && signal_miss==0 && only_ang_miss==0
        est=Fn_TargetPoTrack(est_R,az_beam+hat_azi,el_beam+hat_elv);
        x_est=est(1); y_est=est(2);
        plot([0 max_range*cos((az_beam+theta3/2)*pi/180)],...
            [0 max_range*sin((az_beam+theta3/2)*pi/180)],'Color','blue');
        plot([0 max_range*cos((az_beam-theta3/2)*pi/180)],...
            [0 max_range*sin((az_beam-theta3/2)*pi/180)],'Color','blue');
        %             if signal_miss==0,
        plot(x_est, y_est, 'bo','MarkerFaceColor','b','MarkerSize',7);
        plot([x_est max_range],[y_est max_range],':','Color','blue')
        text(max_range-10,max_range+5,'estimation','Color','blue');
        %             end
        plot(target(k,1),target(k,2),'rx','LineWidth',2);
    end
%     text(-max_range/2,-max_range-10,'Tracking','Color','black','FontSize',14);

subplot(232) %%%%%%%%%%%%
    Fn_PlotMonopulse(theta_s,theta3,az_beam,el_beam);
    if rflag==0 && signal_miss==0 && only_ang_miss==0
        plot(az_beam+[azi_angle -theta3*0.7], el_beam+[elv_angle theta3*0.7], ...
            ':','Color','red')
        text(az_beam-theta3*0.8, el_beam+theta3*0.8,'target','Color','red');
        plot(az_beam+hat_azi, el_beam+hat_elv, 'bo','MarkerFace','b','MarkerSize',10),
        plot(az_beam+[hat_azi theta3*0.7], el_beam+[hat_elv theta3*0.7], ...
            ':','Color','blue')
        text(az_beam+theta3*0.6, el_beam+theta3*0.8,'estimation','Color','blue');
        plot(azim, elev, 'rx','LineWidth',2);
        str=['beam center: (',num2str(round(az_beam*10)/10),...
            ', ',num2str(round(el_beam*10)/10),')'];
        text(az_beam-theta3*0.8, el_beam-theta3*0.8,str,'Color','black');
        az_beam=az_beam+hat_azi;
        el_beam=el_beam+hat_elv;
    elseif signal_miss==1 || only_ang_miss==1
        text(az_beam+-theta3+3,el_beam,'target missing','Color','red','FontSize',14);
    end

graph_signal=gc_signal(1:Nrcvp);
spec_signal=gc_signal;
rmm=1.2;
% rmm=max(abs(graph_signal))*1.2; 
subplot(234),
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
        if (r_gate(1)+Ntx-1)>Nrcvp
            plot(t_rcv, [zeros(1,r_gate(1)-1),ones(1,Nrcvp-r_gate(1)+1)]);
        else
            plot(t_rcv, [zeros(1,r_gate(1)-1),ones(1,Ntx), ...
                zeros(1,Nrcvp-r_gate(1)-Ntx+1)]);
        end
    end
subplot(235),    
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
    text(BW*2, 10^(0), ['v = ',num2str(vel),' m/s']);
    text(BW*2, 10^(-1)*5, ['v_{est} = ',num2str(est_vel),' m/s']);

%% *** Search/Tracking Results ***
if rflag==1
    str_est_R='추적 범위 밖'; str_err_R='';
    str_est_az='추적 범위 밖'; str_err_az='';
    str_est_el='추적 범위 밖'; str_err_el='';
    str_est_V='추적 범위 밖'; str_err_V='';  
elseif signal_miss==1
    str_est_R='Signal Missing'; str_err_R='';
    str_est_az='Signal Missing'; str_err_az='';
    str_est_el='Signal Missing'; str_err_el='';
    str_est_V='Signal Missing'; str_err_V='';
elseif only_ang_miss==1
    str_est_R='Signal Missing'; str_err_R='';
    str_est_az='Signal Missing'; str_err_az='';
    str_est_el='Signal Missing'; str_err_el='';
    str_est_V='Signal Missing'; str_err_V='';
else
    est_az=az_beam+hat_azi;
    est_el=el_beam+hat_elv;
    [R,est_R ; azim,est_az ; elev,est_el ; vel,est_vel]
    R_sqerr_sum=R_sqerr_sum+(R-est_R)^2;
    Az_sqerr_sum=Az_sqerr_sum+(azi_angle-hat_azi)^2;
    Ev_sqerr_sum=Ev_sqerr_sum+(elv_angle-hat_elv)^2;
    V_sqerr_sum=V_sqerr_sum+(vel-est_vel)^2;
    str_est_R=[num2str(est_R),' km']; str_err_R=[num2str(R-est_R),' km'];
    str_est_az=[num2str(est_az),' deg']; str_err_az=[num2str(azi_angle-hat_azi),' deg'];
    str_est_el=[num2str(est_el),' deg']; str_err_el=[num2str(elv_angle-hat_elv),' deg'];
    str_est_V=[num2str(est_vel),' m/s']; str_err_V=[num2str(vel-est_vel),' m/s'];   
end
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['추정: ',str_est_R],'FontSize',10,...
               'Position',[.1 .77 .8 .05],'HorizontalAlignment','left'); % 82
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['오차: ',str_err_R],'FontSize',10,...
               'Position',[.1 .72 .8 .05],'HorizontalAlignment','left');
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['추정: ',str_est_az],'FontSize',10,...
               'Position',[.1 .55 .8 .05],'HorizontalAlignment','left'); % 60
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['오차: ',str_err_az],'FontSize',10,...
               'Position',[.1 .50 .8 .05],'HorizontalAlignment','left');           
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['추정: ',str_est_el],'FontSize',10,...
               'Position',[.1 .33 .8 .05],'HorizontalAlignment','left'); % 38
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['오차: ',str_err_el],'FontSize',10,...
               'Position',[.1 .28 .8 .05],'HorizontalAlignment','left');
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['추정: ',str_est_V],'FontSize',10,...
               'Position',[.1 .11 .8 .05],'HorizontalAlignment','left'); % 16
uicontrol(p2,'Style','text','Units','normalized',...
               'String',['오차: ',str_err_V],'FontSize',10,...
               'Position',[.1 .06 .8 .05],'HorizontalAlignment','left');


end %% ccc

pause(t_intv/proc_speed);
end %% end of num
R_rms=sqrt(R_sqerr_sum/num)     % rms error of estimated range [km]
Az_rms=sqrt(Az_sqerr_sum/num)   % rms error of estimated azim [deg]
Ev_rms=sqrt(Ev_sqerr_sum/num)   % rms error of estimated elev [deg]
V_rms=sqrt(V_sqerr_sum/num)     % rms error of estimated vel [m/s]
p3 = uipanel('Position',[.66 .03 .25 .2]);
uicontrol(p3,'Style','text','Units','normalized',...
               'String','<시뮬레이션 결과>','FontSize',12,...
               'Position',[.12 .7 .7 .2]);
uicontrol(p3,'Style','text','Units','normalized',...
               'String',['rms 거리 오차: ',num2str(R_rms),' km'],'FontSize',10,...
               'Position',[.1 .5 .7 .15],'HorizontalAlignment','left');   
uicontrol(p3,'Style','text','Units','normalized',...
               'String',['rms 방위각 오차: ',num2str(Az_rms),' deg'],'FontSize',10,...
               'Position',[.1 .35 .7 .15],'HorizontalAlignment','left');  
uicontrol(p3,'Style','text','Units','normalized',...
               'String',['rms 고도각 오차: ',num2str(Ev_rms),' deg'],'FontSize',10,...
               'Position',[.1 .2 .7 .15],'HorizontalAlignment','left');   
uicontrol(p3,'Style','text','Units','normalized',...
               'String',['rms 속도 오차: ',num2str(V_rms),' m/s'],'FontSize',10,...
               'Position',[.1 .05 .7 .15],'HorizontalAlignment','left');  