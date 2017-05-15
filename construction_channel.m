

%This code is Implementong the MM wave channel.No output is there in this
%code



clear;clc;
%------------------------System Parameters---------------------------------
Num_BS_Antennas=64; % BS antennas
BSAntennas_Index=0:1:Num_BS_Antennas-1; % Indices of the BS Antennas
Num_BS_RFchains=16; % BS RF chains

Num_MS_Antennas=32; % MS antennas
MSAntennas_Index=0:1:Num_MS_Antennas-1; % Indices of the MS Antennas
Num_MS_RFchains=8;  % MS RF chains

Num_Qbits=7;  % Number of phase shifters quantization bits

% % ---------------------Channel Parameters ---------------------------------




Num_paths=3; % Number of channel paths

%path-loss calculation
Carrier_Freq=28*10^9; % Carrier frequency
lambda=3*10^8/Carrier_Freq; % Wavelength
n_pathloss=3; % Pathloss constant
Tx_Rx_dist=50; % Distance between BS and MS
ro=((lambda/(4*pi*5))^2)*(5/Tx_Rx_dist)^n_pathloss; % Pathloss
Pt_avg=10^(.7); % Average total transmitted power
Pr_avg=Pt_avg*ro; % Average received power

% Noise calculations
Bandwidth=100*10^6;   % Channel bandwidth
No_dB=-173+10*log10(Bandwidth); % Noise power in dB
No=10^(.1*No_dB);    % Noise power absolute




%---------------- Channel Estimation Algorithm Parameters------------------


G_BS=192; % Required resolution for BS AoD
G_MS=192; % Required resolution for MS AoA

K_BS=2;  % Number of Beamforming vectors per stage 
K_MS=2;  % Number of Measurement vectors per stage

Num_paths_est=Num_paths; % Number of desired estimated paths  

S=floor(max(log(G_BS/Num_paths_est)/log(K_BS),log(G_MS/Num_paths_est)/log(K_MS))); % Number of iterations

% Optimized power allocation
Pr_alloc=power_allocation(Num_BS_Antennas,Num_BS_RFchains,BSAntennas_Index,G_BS,G_MS,K_BS,K_MS,Num_paths_est,Num_Qbits);
Pr=abs(Pr_avg*Pr_alloc*S);

%---------------------- Simulation Parameters-------------------------------

ITER=20; % Number of independent realizations (to be averaged over)
SNR_dBa=-40:5:0; % SNR range in dB (for beamforming - after channel estimation)
% Beamsteering vectors generation
for g=1:1:G_BS
    AbG(:,g)=sqrt(1/Num_BS_Antennas)*exp(1j*(2*pi)*BSAntennas_Index*((g-1)/G_BS));
end

% Am generation
for g=1:1:G_MS
    AmG(:,g)=sqrt(1/Num_MS_Antennas)*exp(1j*(2*pi)*MSAntennas_Index*((g-1)/G_MS));
end

%--------------------------------------------------------------------------

for iter=1:1:ITER
    iter 
% Channel parameters (angles of arrival and departure and path gains)
AoD=2*pi*rand(1,Num_paths);
AoA=2*pi*rand(1,Num_paths);
alpha=(sqrt(1/2)*sqrt(1/Num_paths)*(randn(1,Num_paths)+1j*randn(1,Num_paths)));



% Channel construction
Channel=zeros(Num_MS_Antennas,Num_BS_Antennas);
for l=1:1:Num_paths
Abh(:,l)=sqrt(1/Num_BS_Antennas)*exp(1j*BSAntennas_Index*AoD(l));
Amh(:,l)=sqrt(1/Num_MS_Antennas)*exp(1j*MSAntennas_Index*AoA(l));
Channel=Channel+sqrt(Num_BS_Antennas*Num_MS_Antennas)*alpha(l)*Amh(:,l)*Abh(:,l)';
end

% Optimal SVD precoders
[U_H S_H V_H]=svd(Channel);
F_opt=V_H(:,1:Num_paths_est);
W_opt=U_H(:,1:Num_paths_est);
end




