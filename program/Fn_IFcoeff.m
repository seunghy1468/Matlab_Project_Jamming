function IF_result=Fn_IFcoeff(wav,amp_op,f_IF,IF_bw,IF_type,IF_perf,PW,Fs)

switch wav,
    case 1, % cos
        switch IF_bw, % 1-wide; 2-narrow
            case 1, Wn=[f_IF-(1/PW)*3 f_IF+(1/PW)*3];                
            case 2, Wn=[f_IF-(1/PW)*1 f_IF+(1/PW)*1];
        end
        IF_bw_min = Wn(1);
        IF_bw_max = Wn(2);        
        switch IF_type, % 1-butterworth; 2-chevyshev
            case 1, BPF='butter';
            case 2, BPF='cheby';
        end
        switch IF_perf, % 1-good; 2-normal
            case 1, Ord=2;
            case 2, Ord=1;
        end
    case 2, % pulse
        switch IF_bw, % 1-wide; 2-narrow
            case 1, Wn=(1/PW)*3;                
            case 2, Wn=(1/PW)*1;
        end
        IF_bw_min = 0;
        IF_bw_max = Wn;        
        switch IF_type, % 1-butterworth; 2-chevyshev
            case 1, LPF='butter';
            case 2, LPF='cheby';
        end
        switch IF_perf, % 1-good; 2-normal
            case 1, Ord=2;
            case 2, Ord=1;
        end
    case 3, % chirp
        switch IF_bw, % 1-wide; 2-narrow
            case 1, Wn=2*f_IF*1.2;                
            case 2, Wn=2*f_IF*0.8;
        end
        IF_bw_min = 0;
        IF_bw_max = Wn;        
        switch IF_type, % 1-butterworth; 2-chevyshev
            case 1, LPF='butter';
            case 2, LPF='cheby';
        end
        switch IF_perf, % 1-good; 2-normal
            case 1, Ord=5;
            case 2, Ord=2;
        end        
end

passripp = 0.5; % R dB of peak-to-peak ripple in the passband
switch wav,
    case 1, % cos        
        switch BPF,
            case 'butter',  [b,a]=butter(Ord, Wn/Fs); % BPF
            case 'cheby',   [b,a]=cheby1(Ord, passripp, Wn/Fs); % BPF
        end
    otherwise, % pulse or chirp
        switch LPF,
            case 'butter',  [b,a]=butter(Ord, Wn/Fs, 'low'); % LPF
            case 'cheby',   [b,a]=cheby1(Ord, passripp, Wn/Fs, 'low'); % LPF
        end
end

IF_result=[length(a),a,length(b),b];