function G=Fn_Reflection(freq, psi_g, ground_factor, rainfall, tmp)
% ground_factor
% 1 : soil
% 2 : see water
% rainfall, mm/hr
% tmp : temperature, deg
t_s_e1=[2.9    6      10.5    16.7 ; ...  % freq<5
        2.8    5.8    10.3    15.3 ; ...  % freq<10
        2.8    5.6     9.4    12.6 ;...   % freq<15
        2.6    4.9     7.7     9.6];      % freq>=15
t_s_e2=[0.027  0.45   0.75    1.2 ; ...  % freq<5
        0.032  0.87    2.5    4.1 ; ...  % freq<10
        0.035   1.14    3.7    6.3 ;...   % freq<15
        0.03   1.15    4.8    8.5];      % freq>=15  
t_w_e1=[77      75.2    72.3 ; ...  % freq<1
        71      72.1    70.5 ; ...  % freq<4
        56.5    63.2    65.4 ; ...  % freq<8
        47      56.2    60.8];      % freq>=8
t_w_e2=[59.4    73.8    90 ;...     % freq<1
        38.4    38.4    40.2 ; ...  % freq<4
        42      39      36 ; ...    % freq<8
        42.8    40.5    36];        % freq>=8
    
switch ground_factor
    case 1
        if freq<5, row=1;
        elseif freq<10, row=2;
        elseif freq<15, row=3;
        else row=4;
        end
        if rainfall==0, col=1;
        elseif rainfall<10, col=2;
        elseif rainfall<25, col=3;
        else col=4;
        end
        eps1=t_s_e1(row,col); 
        eps2=t_s_e2(row,col);
    case 2
        if freq<1, row=1;
        elseif freq<4, row=2;
        elseif freq<8, row=3;
        else row=4;
        end
        if tmp<5, col=1;
        elseif tmp<15, col=2;
        else col=3;
        end
        eps1=t_w_e1(row,col); 
        eps2=t_w_e2(row,col);
end
j=sqrt(-1);
e=eps1-j*eps2;
A1=sin(psi_g);
A2=sqrt(e-(cos(psi_g)).^2);
Gv=abs((e*A1-A2)./(e*A1+A2));
Gh=abs((A1-A2)./(A1+A2));
G=[Gv,Gh];
end