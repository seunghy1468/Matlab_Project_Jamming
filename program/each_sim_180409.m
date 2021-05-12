clear all; 
clf
pause(1)

%%--- 전자전 위협체
type_radar=2; % 1-search, 2-monopulse, 3-monopulse_doppler, 4-CW_doppler
% if type_radar==1 (search), wav=1,2 (rect pulse), PRI=400
% if type_radar==2 (only monopulse), wav=1,2,or3 (all), PRI=200
% if type_radar==3 (monopulse + pulse doppler), wav=1 (cos pulse), PRI=200
% if type_radar==4 (monopulse + CW doppler), wav=1 (cos wave), PRI=PW
% if type_radar==5 (search CW doppler), wav=1 (cos wave), PRI=PW

%%--- 재머
type_jam=41; 
% 0 -no jamming
% 11-barrage_noise, 12-spot_noise
% 21-RGPO, 22-RGPI
% 31-VGPO, 32-VGPI
% 41-plus_angle, 42-minus_angle
% 61-V/RGPI, 62-V/RGPO

%% ===== GUI inputs =====

switch type_radar
    case 1 %%% search
        
        %%--- 수집변수
        wav=1; % 1-cos, 2-pulse, 3-chirp
        PW=50; % [micro sec] pulse width
        PP=70; % [dBm], transmit pulse power
        PRI=400; % [micro seq] pulse repetition interval : 400 -> max range 60km
        f_RF=1*1000; % [MHz] carrier frequency
        Tscan=12; % [s] scan period
        
        %%--- 모의실험 환경
        proc_time=76; % [s] processing time : search->200, monopulse->20
        t_intv=1; % [s] processing interval
        num=proc_time/t_intv; % number of simulation > stt_jam+dur_jam (=12)
        Fs=4000/PRI; % [Msps] sampling rate, 4000 samples/PRI
        proc_speed=5; % 1,2,5,10
        
        %%--- 전자전 수신기 (항공기)
        RCS=20; % [m^2] radar cross section
        vel=300; % [m/s] velocity
        ip_target=[-44 43]; % [km] initial position of target : search_[-40 50], track_[-23 20]
        ad_target=350; % [deg] angle of direction : search_355, track_330
        
        %%-- 전자전 위협체
        Pfa = 10^(-5); % false alarm probability
        Pd = 0.99999; % detection probability
        FBank=2^16; % number of filter = fft size
        theta_s=3; % [deg] squint angle
        %%-- H/W (Antenna)
        theta3=360/Tscan; % [deg] 3dB beamwidth
        Gt=40; % [dB] tx antenna gain
        Gr=40; % [dB] rx antenna gain
        %%-- H/W (AGC / Limiter)
        gain_control=1; % 1-Limiter, 2-AGC
        agc_dyn=1; % 1-wide; 2-narrow
        agc_delay=4; % [s] 안정화까지 걸리는 시간
        %%- H/W (IF, Mixer)
        f_MF=f_RF/1000; % [MHz] middle frequency
        f_IF=500; % [kHz] intermediate frequency : 2/PW (2 periods/PW) <= f_IF << Fs/2
        G_IF=30; % [dB]
        amp_op=1; % 1-wide; 2-narrow
        IF_bw=1; % 1-wide; 2-narrow
        IF_type=1; % 1-butterworth; 2-chevyshev
        IF_perf=1; % 1-good; 2-normal
        switch amp_op
            case 1, op_range=30; % [dB]
            case 2, op_range=20; % [dB]
        end
        
    case 2 %%% monopulse
        
        %%--- 수집변수
        wav=1; % 1-cos, 2-pulse, 3-chirp
        PW=20; % [micro sec] pulse width
        PP=70; % [dBm], transmit pulse power
        PRI=200; % [micro seq] pulse repetition interval : 400 -> max range 60km
        f_RF=1*1000; % [MHz] carrier frequency
        Tscan=0; % [s] scan period
        
        %%--- 모의실험 환경
        proc_time=30; % [s] processing time : search->200, monopulse->20
        t_intv=0.5; % [s] processing interval
        num=proc_time/t_intv; % number of simulation > stt_jam+dur_jam (=12)
        Fs=4000/PRI; % [Msps] sampling rate, 4000 samples/PRI
        proc_speed=5; % 1,2,5,10
        
        %%--- 전자전 수신기 (항공기)
        RCS=20; % [m^2] radar cross section
        vel=300; % [m/s] velocity
        ip_target=[-15 15]; % [km] initial position of target : search_[-40 50], track_[-23 20]
        ad_target=355; % [deg] angle of direction : search_355, track_330
        
        %%-- 전자전 위협체
        Pfa = 10^(-5); % false alarm probability
        Pd = 0.99999; % detection probability
        FBank=2^16; % number of filter = fft size
        theta_s=3; % [deg] squint angle
        %%-- H/W (Antenna)
        theta3=10; % [deg] 3dB beamwidth
        Gt=40; % [dB] tx antenna gain
        Gr=40; % [dB] rx antenna gain
        %%-- H/W (AGC / Limiter)
        gain_control=2; % 1-Limiter, 2-AGC
        agc_dyn=1; % 1-wide; 2-narrow
        agc_delay=4; % [s] 안정화까지 걸리는 시간
        %%- H/W (IF, Mixer)
        f_MF=f_RF/1000; % [MHz] middle frequency
        f_IF=500; % [kHz] intermediate frequency : 2/PW (2 periods/PW) <= f_IF << Fs/2
        G_IF=30; % [dB]
        amp_op=1; % 1-wide; 2-narrow
        IF_bw=1; % 1-wide; 2-narrow
        IF_type=1; % 1-butterworth; 2-chevyshev
        IF_perf=1; % 1-good; 2-normal
        switch amp_op
            case 1, op_range=30; % [dB]
            case 2, op_range=20; % [dB]
        end
        
    case 3 %%% monopulse doppler
        
        %%--- 수집변수
        wav=1; % 1-cos, 2-pulse, 3-chirp
        PW=50; % [micro sec] pulse width
        PP=70; % [dBm], transmit pulse power
        PRI=200; % [micro seq] pulse repetition interval : 400 -> max range 60km
        f_RF=1*1000; % [MHz] carrier frequency
        Tscan=0; % [s] scan period
        
        %%--- 모의실험 환경
        proc_time=20; % [s] processing time : search->200, monopulse->20
        t_intv=0.5; % [s] processing interval
        num=proc_time/t_intv; % number of simulation > stt_jam+dur_jam (=12)
        Fs=4000/PRI; % [Msps] sampling rate, 4000 samples/PRI
        proc_speed=5; % 1,2,5,10
        
        %%--- 전자전 수신기 (항공기)
        RCS=20; % [m^2] radar cross section
        vel=300; % [m/s] velocity
        ip_target=[-15 15]; % [km] initial position of target : search_[-40 50], track_[-23 20]
        ad_target=350; % [deg] angle of direction : search_355, track_330
        
        %%-- 전자전 위협체
        Pfa = 10^(-5); % false alarm probability
        Pd = 0.99999; % detection probability
        FBank=2^16; % number of filter = fft size
        theta_s=3; % [deg] squint angle
        %%-- H/W (Antenna)
        theta3=10; % [deg] 3dB beamwidth
        Gt=40; % [dB] tx antenna gain
        Gr=40; % [dB] rx antenna gain
        %%-- H/W (AGC / Limiter)
        gain_control=2; % 1-Limiter, 2-AGC
        agc_dyn=1; % 1-wide; 2-narrow
        agc_delay=4; % [s] 안정화까지 걸리는 시간
        %%- H/W (IF, Mixer)
        f_MF=f_RF/1000; % [MHz] middle frequency
        f_IF=500; % [kHz] intermediate frequency : 2/PW (2 periods/PW) <= f_IF << Fs/2
        G_IF=30; % [dB]
        amp_op=1; % 1-wide; 2-narrow
        IF_bw=1; % 1-wide; 2-narrow
        IF_type=1; % 1-butterworth; 2-chevyshev
        IF_perf=1; % 1-good; 2-normal
        switch amp_op
            case 1, op_range=30; % [dB]
            case 2, op_range=20; % [dB]
        end
        
    case 4 %%% Monopulse CW doppler
        
        %%--- 수집변수
        wav=1; % 1-cos, 2-pulse, 3-chirp
        PW=200; % [micro sec] pulse width
        PP=70; % [dBm], transmit pulse power
        PRI=200; % [micro seq] pulse repetition interval : 400 -> max range 60km
        f_RF=1*1000; % [MHz] carrier frequency
        Tscan=0; % [s] scan period
        
        %%--- 모의실험 환경
        proc_time=20; % [s] processing time : search->200, monopulse->20
        t_intv=0.5; % [s] processing interval
        num=proc_time/t_intv; % number of simulation > stt_jam+dur_jam (=12)
        Fs=4000/PRI; % [Msps] sampling rate, 4000 samples/PRI
        proc_speed=5; % 1,2,5,10
        
        %%--- 전자전 수신기 (항공기)
        RCS=20; % [m^2] radar cross section
        vel=600; % [m/s] velocity
        ip_target=[-15 15]; % [km] initial position of target : search_[-40 50], track_[-23 20]
        ad_target=350; % [deg] angle of direction : search_355, track_330
        
        %%-- 전자전 위협체
        Pfa = 10^(-5); % false alarm probability
        Pd = 0.99999; % detection probability
        FBank=2^16; % number of filter = fft size
        theta_s=3; % [deg] squint angle
        %%-- H/W (Antenna)
        theta3=10; % [deg] 3dB beamwidth
        Gt=40; % [dB] tx antenna gain
        Gr=40; % [dB] rx antenna gain
        %%-- H/W (AGC / Limiter)
        gain_control=1; % 1-Limiter, 2-AGC
        agc_dyn=1; % 1-wide; 2-narrow
        agc_delay=4; % [s] 안정화까지 걸리는 시간
        %%- H/W (IF, Mixer)
        f_MF=f_RF/1000; % [MHz] middle frequency
        f_IF=500; % [kHz] intermediate frequency : 2/PW (2 periods/PW) <= f_IF << Fs/2
        G_IF=30; % [dB]
        amp_op=1; % 1-wide; 2-narrow
        IF_bw=1; % 1-wide; 2-narrow
        IF_type=1; % 1-butterworth; 2-chevyshev
        IF_perf=1; % 1-good; 2-normal
        switch amp_op
            case 1, op_range=30; % [dB]
            case 2, op_range=20; % [dB]
        end
        
    case 5 %%% Search CW doppler
        
        %%--- 수집변수
        wav=1; % 1-cos, 2-pulse, 3-chirp
        PW=400; % [micro sec] pulse width
        PP=70; % [dBm], transmit pulse power
        PRI=400; % [micro seq] pulse repetition interval : 400 -> max range 60km
        f_RF=1*1000; % [MHz] carrier frequency
        Tscan=12; % [s] scan period
        
        %%--- 모의실험 환경
        proc_time=5; % [s] processing time : search->200, monopulse->20
        t_intv=1; % [s] processing interval
        num=proc_time/t_intv; % number of simulation > stt_jam+dur_jam (=12)
        Fs=4000/PRI; % [Msps] sampling rate, 4000 samples/PRI
        proc_speed=5; % 1,2,5,10
        
        %%--- 전자전 수신기 (항공기)
        RCS=20; % [m^2] radar cross section
        vel=300; % [m/s] velocity
        ip_target=[-20 20]; % [km] initial position of target : search_[-40 50], track_[-23 20]
        ad_target=350; % [deg] angle of direction : search_355, track_330
        
        %%-- 전자전 위협체
        Pfa = 10^(-5); % false alarm probability
        Pd = 0.99999; % detection probability
        FBank=2^16; % number of filter = fft size
        theta_s=3; % [deg] squint angle
        %%-- H/W (Antenna)
        theta3=360/Tscan; % [deg] 3dB beamwidth
        Gt=40; % [dB] tx antenna gain
        Gr=40; % [dB] rx antenna gain
        %%-- H/W (AGC / Limiter)
        gain_control=1; % 1-Limiter, 2-AGC
        agc_dyn=1; % 1-wide; 2-narrow
        agc_delay=4; % [s] 안정화까지 걸리는 시간
        %%- H/W (IF, Mixer)
        f_MF=f_RF/1000; % [MHz] middle frequency
        f_IF=500; % [kHz] intermediate frequency : 2/PW (2 periods/PW) <= f_IF << Fs/2
        G_IF=30; % [dB]
        amp_op=1; % 1-wide; 2-narrow
        IF_bw=1; % 1-wide; 2-narrow
        IF_type=1; % 1-butterworth; 2-chevyshev
        IF_perf=1; % 1-good; 2-normal
        switch amp_op
            case 1, op_range=30; % [dB]
            case 2, op_range=20; % [dB]
        end
end

%%--- 재머
BT_range=10; % [km] burn through range (J=S)
Gj=30; % [dB] jammer antenna gain
dur_jamming=30; % [s] duration of jamming
dur_jam=dur_jamming/t_intv;
delay_jamming=5; % [s] detected -> wait during delay_jam -> start jamming
switch type_jam
    case 11 % barrrage noise
        wc1=0; wc2=3; % noise band [MHz]
    case 12 % spot noise
        spot_num=5;
    case 21
        vel_jam=300; % [m/s] = velocity
    case 22
        vel_jam=300; % [m/s] = velocity
    case 31 %% VGPO 
        acc_jam=20; % [m/(s^2)] = acceleration
    case 32 %% VGPI
        acc_jam=20; % [m/(s^2)]
    case 61 %% VGPO+RGPO
        acc_jam=20; % [m/(s^2)] = acceleration
    case 62 %% VGPI+RGPI
        acc_jam=20; % [m/(s^2)]        
end

%%--- 전자전 송수신환경
temper=10; % [deg] temperature
ground_type=1; % 1-soil(육지), 2-sea(바다)
rainfall=0; % [mm/hr] 

%%--- 상수 (미입력)
c=3; % [*10^8 m/s] velocity of light
h_radar=0.3; % [km] 레이더 높이 
h_target=3; % [km] 수신기 높이 (z 좌표)
g0=1; % antenna center gain

save each_sim_2018.mat

switch type_radar
    case 1
        Fn_RadarSearchPulse();
    case 2
        Fn_RadarMonoPulse();        
    case 3
        Fn_RadarMonoPulseDoppler();
    case 4
        Fn_RadarMonoContWave();
    case 5
        Fn_RadarSearchContWave();
end