function [E, E_bit, X, tau, Reff, Pirr_dBm] = comm_consumption(pc,ps,pg)
% ----------------------------------------------------------------------
%  Function which calculates the energy consumed for sensing audio data. 
%  Author: Fernando Rosas and Gert Dekkers, Imperial College London/KU Leuven
% ----------------------------------------------------------------------
% Syntaxis: comm_consumption(pc,ps,pg)
% Inputs:
% (1,2,3) pc,ps,pg      struct containing all params of communication, sensing and general resp.
% Outputs:
% (1) E - total energy consumed by the transceiver/receiver [Etx,Erx] [mJ]
% (2) E_bit - total energy consumed by the transceiver/receiver per data bit [Etx,Erx] [mJ]
% (3) X - vector which contains the percentage of the total consumption of
% the transceiver corresponding to X(1): startup, X(2): encoding, X(2): electronic components involved
% in the forward transmisison, X(3): electronic components involved in the
% feedback transmissions, X(4): PA.
% (4) Tau - average number of transmission trials until a frame is
% correctly decoded in the receiver
% (5) Reff - Effective data rate [bits/s]    
% (6) Pirr_dBm - Irradiated power [dBm]

% load parameters
v2struct(pc); % expand field to workspace
v2struct(ps); % expand field to workspace
v2struct(pg); % expand field to workspace

% general things
b = log2(M);                                                                % number of bits per simbol
omega = 1;                                                                  % no multiplexion


%% Error correcting code
if ~Ladj
    % get matching BCD param
    [ntx, Ktx,ttx] = BCDparam(Ltx,ttx,N_T); 
    [nrx, Krx,trx] = BCDparam(Lrx,trx,N_R);
else
    L_list = [7 15 31 63 127 255 511];
    for k = 1:length(L_list) % get an adjusted payload size that better mathes N_T
        [ntx_tmp, Ktx_tmp,ttx_tmp] = BCDparam(L_list(k),ttx,N_T); % get matching BCD param
        if ~isempty(ttx_tmp) && N_T<=Ktx_tmp, break; end; 
    end
    ntx = ntx_tmp; Ktx = Ktx_tmp; ttx = ttx_tmp; Ltx = ntx+1;
    for k = 1:length(L_list)
        [nrx_tmp, Krx_tmp,trx_tmp] = BCDparam(L_list(k),trx,N_R); % get matching BCD param
        if ~isempty(trx_tmp) && N_R<=Krx_tmp, break; end; 
    end
    nrx = nrx_tmp; Krx = Krx_tmp; trx = trx_tmp; Lrx = nrx+1;
end
% adjust usefull bytes to integer frame multiple and integer byte multple
N_T = ceil(N_T/Ktx)*Ktx; if isnan(N_T), N_T=0; end;
N_R = ceil(N_R/Krx)*Krx; if isnan(N_R), N_R=0; end;
% get coderate
rtx = Ktx/ntx;
rrx = Krx/nrx;
% consumption
Eenc = 0;
% (for now we just neglect the Tx encoding for BCH codes) - section 5.3
Nadd  = (2*trx-1)*trx + 2*trx^2;
Nprod =     2*trx*trx + 2*trx^2;
Edec =  Eop * (Nadd * c(addid) + Nprod * c(multid));

%% times per bit
Tb_tx = 1/(rtx*Rs) * (1/b + H/Ltx + (Oa+Ob)/Ltx);                               % [s] - eq. 17
Tfb_tx = F/(rtx*Rs*Ltx);                                                        % [s] - eq. 18
Tb_rx = 1/(rrx*Rs) * (1/b + H/Lrx + (Oa+Ob)/Lrx);                               % [s] - eq. 17
Tfb_rx = F/(rrx*Rs*Lrx);                                                        % [s] - eq. 18
if isinf(Tb_tx), Tb_tx = 0; end;
if isinf(Tfb_tx), Tfb_tx = 0; end;
if isinf(Tb_rx), Tb_rx = 0; end;
if isinf(Tfb_rx), Tfb_rx = 0; end;

%% consumption of electric components
Pdac_tx = .5* ( Vdd_dac*I0_dac*(2^n1 -1) + n1*Cp_dac*fs_DAC*Vdd_dac^2);     % [mW] - eq. 34
Padc_rx =  2^n2 * fs_ADC * FOMadc;                                          % [mW] - eq. 10
Petx = Pdac_tx + 2*Pfilter + Plo + Pmixer;                                  % [mW] - eq. 31 (partial)
Perx  = 3*Pfilter + Plna_rx + Plo + Pmixer + Pvga + Padc_rx;                % [mW] - eq. 32 (partial)
Eetx_b = Petx * Tb_tx;                                                         % [mJ]
Eetx_fb = Petx * Tfb_tx;                                                       % [mJ]
Eerx_b = Perx * Tb_rx;                                                         % [mJ]
Eerx_fb = Perx * Tfb_rx;                                                       % [mJ]

%% Calculating the constant "A" and PA consumption
if M == 2                                                                   % Peak-to-average power ratio
    PAPR=1;
else
    PAPR =  3 * (sqrt(M)-1)/(sqrt(M)+1);
end

c=3*10^8;                                                                   % Speed of light [m/s]
lambda=c/fc;                                                                % Wavelength [m]
A0 = 1/(Gt*Gr) * (4*pi/lambda);                                             % Free space attenuation at 1 meter - eq. 25

if PA == 0     % if ClassA PA - eq. 21
    beta = 1;
    eta_eff  = 0.5;
elseif PA == 1 % if ClassB PA - eq. 22
    beta = 0.5;
    eta_eff  = 0.785;
end

k = 1.3806488 * 10^-23;                                                     % Boltzmann constant [m^2 kg / s^2 K]
T_K = T_C+273;                                                              % Room temperature [K]
N0 = 10*log10(k*T_K*10^3);                                                  % Noise spectral density [dBm]
Pnoise = N0 + 10*log10(W) + Nf + Ml;                                        % Total noise power [dBm] - section 5.1
Pnoise_abs = 10^(Pnoise/10);                                                % Total noise power [mW] - section 5.1
A = (PAPR/Spa)^beta * Pnoise_abs * A0 / eta_eff;                            % [mW]

EbN0_abs = 10^(e/10); 
SNR_abs = b*EbN0_abs;
Ppa = A * omega * d^a * SNR_abs;                                            % [mW]
Epa_tx_b = Ppa * Tb_tx;                                                       	% [mJ]
Epa_rx_fb = Ppa * Tfb_rx;                                                      % [mJ]

%% Retransmission statistics
ws = v2struct;
[qx_tx, tau_tx] = errors(ws,ntx,ttx); % section 5.4
[qx_rx, tau_rx] = errors(ws,nrx,trx); % section 5.4

%% Prepare output
% energy consumption forward information frame
Etx_bit(1) = Est/((1-qx_tx)*N_T) + Eenc + (Eetx_b + Epa_tx_b + Eetx_fb)*tau_tx;      % [mJ] - eq. 14
Etx = Etx_bit*N_T;
% energy consumption backward information frame
Erx_bit = Est/((1-qx_rx)*N_R) + (Edec + Eerx_b + Eerx_fb + Epa_rx_fb)*tau_rx;      % [mJ] - eq. 15
Erx = Erx_bit*N_R;
% effective output bandwidth
Ttotal = (Tb_tx + Tfb_tx)*tau_tx;
Reff = 1/Ttotal; 
% power distribution components forward information frame
Xtx = [Est/((1-qx_tx)*N_T), Eenc,      Eetx_b*tau_tx, Eerx_fb*tau_tx, Epa_tx_b*tau_tx]/Etx_bit;
Xrx = [Est/((1-qx_rx)*N_R), Edec*tau_rx,  Eerx_b*tau_rx, Eetx_fb*tau_rx, Epa_rx_fb*tau_rx]/Erx_bit;
% Irradiated power   
Pirr = (Spa/PAPR)^beta * eta_eff * Ppa;
Pirr_dBm = 10*log10(Pirr);
% Combine
if N_R==0, Erx=0; Erx_bit=0; Xrx = zeros(size(Xrx)); end;
if N_T==0, Etx=0; Etx_bit=0; Xtx = zeros(size(Xtx)); end;
E = [Etx Erx];
E_bit = [Etx_bit Erx_bit];
X = [Xtx; Xrx];
end


function [ qx, tau] = errors(ws,n,t)
v2struct(ws) % expand field
if strcmp(ch_correlation ,'ff')
    mod = 0;
    M_header = 2; % header symbols use BPSK
    Ph = Pbit(e,M_header,channel,K);
    Pb  = Pbit(e,M,channel,K);
        
    P_BCH = 0;
    for j=0:t
        %  Pb_BCH(i,:)=Pb_BCH(i,:)+(factorial(ntx)/(factorial(j)*factorial(ntx-j)))*Pb.^j.*(1-Pb).^(ntx-j);
        P_BCH = P_BCH + nchoosek(n,j) * (1-Pb)^(n-j) * Pb^j;
    end
    
    P_correct_frame = (1 - Ph)^H * P_BCH;
    tau = 1 / P_correct_frame;
    qx = 0; %outage probability   
elseif strcmp(ch_correlation ,'static')
    tau = 1;
    qx = 0; %outage probability
end
end

function P = Pbit(e,M,channel,K)
%
% ---------------------------------------------------------------------
%  Function to calculate the EXACT bit error rate of M-QAM modulations
%
% Pawgn(e,M, channel, K) 
% e: EbN0,
% M: M-ary number,
% channel: 0 -AWGN, 1 - Rayleigh, 3-Nakagami-M,
% K: Nakagami-M parameters (works fine until K larger than 35). Note that
% the performance is much faster if K is an integer.
%

EbN0 = 10^(e/10);

if M == 2
    
    if channel == 0
        P=.5*erfc(sqrt(abs(EbN0)))';
    elseif channel == 1
        P=0.5*(1-sqrt((EbN0)./(1+(EbN0))));

    elseif channel == 3
        
        paso=0.01;
        theta=0:paso:pi/2;
        SNR=EbN0'*ones(size(theta));
        Theta=ones(size(e))'*theta;
        a = 1;
        c = 1;        
        H=K;
        b = log2(M);
        F = ( 1 + a*b*SNR./(H*sin(Theta).^2) ) .^-H;
        P = (c/pi)* sum(F,2) * paso;
    
    end
    
    
else
    
    upbound1 = .5*log2(M);
    Pcases = zeros(upbound1,1);
    
    for k = 1: upbound1
        
        
        upbound2 = ( 1 - 2^-k )* sqrt(M) - 1;
        F = zeros(upbound2+1,1);
        for i = 0: upbound2
            
            arg = i*2^(k-1)/sqrt(M);
            w = (-1)^floor(arg) * ( 2^(k-1) - floor(arg + .5) );
            
            if channel == 0
                
                arg1 = sqrt( 3*log2(M)*EbN0 / ( 2*(M-1) ));
                F(i+1) = w * erfc( (2*i + 1) * arg1);
                
            elseif channel == 1
                
                %arg1 = 3*(2*i+1)^2 *log2(M)*EbN0 / (2*(M-1));
                arg1 = 2*(M-1) / ( 3*(2*i+1)^2 *log2(M)*EbN0 );
                G = sqrt( 1/ (1 + arg1) );
                F(i+1) = w * ( 1 - G);
                
            elseif channel == 3
                
                a = 3 *(2*i+1)^2 * log2(M)/(2*(M-1));
                
                if K == round(K)
                    
                    H = zeros(K,1);
                    for j = 0:K-1
                        arg = K/(a*EbN0);
                        arg1 = (1 + arg)^(-j-.5) +1;
                        H(j+1) = (-1)^j * nchoosek(K-1,j) * arg1 / (2*j+1);
                    end
                    arg0 = sum(H);
                    
                    G = 1 - gamma(K+.5)/(gamma(K)*sqrt(pi)) * arg0;
                    
                else
                    
                    arg0 = a*EbN0/K;
                    arg1 = hypergeom( [.5, K+.5], 3/2, -arg0);
                    G = .5 - gamma(K+.5)/(gamma(K)*sqrt(pi)) * sqrt(arg0) * arg1;
                    
                end
                
                F(i+1) = 2 * w * G;
            end
        end
        Pcases(k) = 1/sqrt(M) * sum(F);
    end
    P = 1/upbound1 * sum(Pcases);
end
end

function [n, K, t] = BCDparam(n,t_sep,N)
    % Based on Lin and Costello - Error Control Coding (Appendix C)
    % get hardcoded values
    if t_sep==0
        n = ceil(N/8)*8;
        K = n;
        t = 0;
    elseif n<=7
        n = 7;
        K = [7 4];
        t = [0 1];
    elseif n<=15
        n = 15;
        K = [15 11 7 5];
        t = [0 1 2 3];
    elseif n<=31
        n = 31;
        K = [31 26 21 16 11 6];
        t = [0 1:3 5 7];
    elseif n<=63
        n = 63;
        K = [63 57 51 45 39 36 30 24 18 16 10 7];
        t = [0 1:7 10 11 13 15];
    elseif n<=127
        n = 127;
        K = [127 120 113 106 99 92 85 78 71 64 57 50 43 36 29 22 15 8];
        t = [0 1:7 9:11 13:15 21 23 27 31];
    elseif n<=255
        n = 255;
        K = [255 247 239 231 223 215 207 199 191 187 179 171 163 155 147 139 131 123 115 107 99 91 87 79 71 63 55 47 45 37 29 21 13 9];
        t = [0 1:15 18 19 21:23 25:27 29:31 42 43 45 47 55 59 63];
    elseif n<=511
        n = 511;
        K = [511 502 493 484 475 466 457 448 439 430 421 412 403 394 385 376 367 358 349 340 331 322 313 304 295 286 277 268 259 250 241 238 229 220 211 202 193 184 175 166 157 148 139 130 121 112 103 94 85 76 67 58 49 40 31 28 19 10];
        t = [0 1:15 17:23 25:31 36:39 41:43 45:47 51 53:55 58 59 61:63 85 87 91 93 95 109 111 119 127];
    elseif n>511
        n = 1023;
        K = [1023,  1013,   1003,   993,    983,    973,    963,    953, 943, 933, 923, 913, 903, 893, 883, 873, 863, 858, 848, 838, 828, 818, 808, 798, 788, 778, 768, 758, 748, 738, 728, 718, 708, 698, 688, 678, 668, 658, 648, 638, 628, 618, 608, 598, 588, 578, 573, 563, 553, 543, 533, 523, 513, 503, 493, 483, 473, 463, 453, 443, 433, 423, 413, 403, 393, 383, 378, 368, 358, 348, 338, 328, 318, 308, 298, 288, 278, 268, 258, 248, 238, 228, 218, 208, 203, 193, 183, 173, 163, 153, 143, 133, 123, 121, 111, 101, 91, 86, 76, 66, 56, 46, 36, 26, 16, 11];
        t = [0  1:31 34:39 41: 47 49:55 57:63 73:75 77:79 82 83 85:87 89:91 93:95 102 103 106:108 110 111 115 117:119 122 123 125 126 127 170 171 173 175 181 183 187 189 191 219 223 239 247 255];
    else
        error('Value of n not supported');
    end
    % select K/t
    ind = find(t>=t_sep); %pick one of atleast t ecc bits
    if ~isempty(ind)
        K = K(ind(1));
        t = t(ind(1));   
    else
        K = 0; t = [];
    end
end