function result=Fn_Waveform(wav,f_IF,fd,PW,tc,rho,flag)
% f_IF [MHz]
% fd [MHz]
% flag==0 : tx
% flag==1 : rx without pulse spreading
% flag==2 : rx with pulse spreading
% f_IF unit == fd unit

switch flag
    case 0
        switch wav
            case 1, result=cos(2*pi*f_IF*tc);
            case 2, result=ones(1,length(tc));
            case 3, result=chirp(tc,0,PW,2*f_IF,'linear');
        end
    case 1
        switch wav
            case 1, result=cos(2*pi*(f_IF+fd)*tc);
            case 2, result=ones(1,length(tc));
            case 3, tr=chirp(tc,0,PW,2*f_IF,'linear');
                result=tr.*cos(2*pi*fd*tc);
        end
    case 2
        tt=tc(2)-tc(1); % [micro sec]
        Ntc=length(tc); % number of samples of tc
        Nsp=(rho/2)/tt; % number of samples in rho/2
        pulse_spread=[1-exp(-10*tc(1:Nsp)/(rho/2)), ...
            ones(1,Ntc-2*Nsp), exp(-10*tc(1:Nsp)/(rho/2))];
        plot(tc,pulse_spread)
        switch wav
            case 1, result=cos(2*pi*(f_IF+fd)*tc).*pulse_spread;
            case 2, result=pulse_spread;
            case 3, tr=chirp(tc,0,PW,2*f_IF,'linear');
                result=tr.*cos(2*pi*fd*tc).*pulse_spread;
        end
end