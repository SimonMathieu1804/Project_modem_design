% LELEC2880: Modem design - Project
% Authors: DE COCK Justin, DELHAYE Quentin, SIMON Mathieu
% Date: 12/04/20

clc;
%%  Step 1 : Basic OFDM chain
%======================= given values =====================================
N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
Nb = 2*N; % block size

%======================= Es=1 and N0=0.1 ==================================
N0 = 0.1;
% 1) vector of 2048 random bits (to send 4 OFDM packets)
bits = randi([0 1],1,4*2*Nb);
% 2) Symbol mapping
map = bits;
map(map==0) = -1;
map = sqrt(2)/2*map;
symbols = zeros(4*Nb,1);
for k=1:4*Nb
    symbols(k)=map(2*k-1)+1i*map(2*k);
end
figure(1);
x = real(symbols); y = imag(symbols);
scatter(x,y,40,'o','filled','r'); title('Tx constellation','Fontsize',16);
xlabel('In phase amplitude','Fontsize',14); ylabel('Quandrature amplitude','Fontsize',14);
% 3) Seriel to parralel
parallel = [symbols(1:Nb) symbols(Nb+1:2*Nb) symbols(2*Nb+1:3*Nb) symbols(3*Nb+1:4*Nb)]; % each column is a block of 256 symbols
% 4) IFFT on the blocks
parallel = ifft(parallel);
% 5) Cyclic prefix insertion
CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4)];
paralel_CP = [CP ; parallel];
% 6) parallel to serial
serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).'];

% Pulse Shapping
alpha=0.2;
N_truncated=10;
u= rcosdesign(alpha,N_truncated,10,'sqrt');
E_u= u*u'; %E_u should be equal to one
fvtool(u,'impulse'); % Plot the filter
fvtool(u,'freq');
x = upfirdn(serial, u, 10);
% Shift at carrier freq ??
% AWGN Channel
x = x + randn(size(x))*2*N0;
% bring back from carrier freq ??
% Matched filter
y = upfirdn(x, u, 1, 10);

% 7) AWGN channel
y = serial+ randn(size(serial))*sqrt(N0)+ randn(size(serial))*sqrt(N0)*1i; % line to comment in order to use pulse shaping
% 8) serial to parralel
y=y.';
parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L))];
% 9) Remove CP
parallelRx = parallelRx((L+1):end,:);
% 10) FFT on the blocks
parallelRx = fft(parallelRx);
% 11) parallel to serial
output = [parallelRx(:,1).' parallelRx(:,2).' parallelRx(:,3).' parallelRx(:,4).'];
figure(4);
x = real(output); y = imag(output);
scatter(x,y,40,'o','filled','r'); title('Rx constellation','Fontsize',16);
xlabel('In phase amplitude','Fontsize',14); ylabel('Quandrature amplitude','Fontsize',14);
% 12) demapping
output_bits = zeros(4*2*Nb,1);
for k=1:4*Nb
    output_bits(2*k-1)=real(output(k));
    output_bits(2*k)=imag(output(k));
end
% 13) decision
output_bits(output_bits<=0)=0;
output_bits(output_bits>0)=1;

BER = sum(output_bits.'~=bits)/2048

%======================= Es=1 and varying N0 ==============================

% BER calculation for each noise level
Nsnr=20;
Es_N0_dB=linspace(0,15,Nsnr);
Es_N0=10.^(Es_N0_dB/10);

BER=zeros(Nsnr,1);
for index_SNR=1:Nsnr
    N0=1/Es_N0(index_SNR);
    %---------------------------
    % 1) vector of 2048 random bits (to send 4 OFDM packets)
    bits = randi([0 1],1,4*2*Nb);
    % 2) Symbol mapping
    map = bits;
    map(map==0) = -1;
    map = sqrt(2)/2*map;
    symbols = zeros(4*Nb,1);
    for k=1:4*Nb
        symbols(k)=map(2*k-1)+1i*map(2*k);
    end
    % 3) Seriel to parralel
    parallel = [symbols(1:Nb) symbols(Nb+1:2*Nb) symbols(2*Nb+1:3*Nb) symbols(3*Nb+1:4*Nb)]; % each column is a block of 256 symbols
    % 4) IFFT on the blocks
    parallel = ifft(parallel);
    % 5) Cyclic prefix insertion
    CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4)];
    paralel_CP = [CP ; parallel];
    % 6) parallel to serial
    serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).'];
    % 7) AWGN channel
    y = serial+ randn(size(serial))*sqrt(N0)+ randn(size(serial))*sqrt(N0)*1i;
    % 8) serial to parralel
    y=y.';
    parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L))];
    % 9) Remove CP
    parallelRx = parallelRx((L+1):end,:);
    % 10) FFT on the blocks
    parallelRx = fft(parallelRx);
    % 11) parallel to serial
    output = [parallelRx(:,1).' parallelRx(:,2).' parallelRx(:,3).' parallelRx(:,4).'];
    % 12) demapping
    output_bits = zeros(4*2*Nb,1);
    for k=1:4*Nb
        output_bits(2*k-1)=real(output(k));
        output_bits(2*k)=imag(output(k));
    end
    % 13) decision
    output_bits(output_bits<=0)=0;
    output_bits(output_bits>0)=1;
    
    %---------------------------
    N_symbols = 4*2*Nb;
    BER(index_SNR)=sum(output_bits.'~=bits)/N_symbols;
end


k=4;
M=2^k;
x=sqrt(3*k*Es_N0/(M-1));
theoretical_BER=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));

figure(5);
semilogy(Es_N0_dB,theoretical_BER,'-r','LineWidth',1.5);
hold on;
semilogy(Es_N0_dB,BER/2,'-xb','LineWidth',1.5,'MarkerSize',8);
grid;
xlabel('E_S/N_0 [dB]'); ylabel('SER'); legend('Theory (4QAM)','Simulated');

%% Step 2: Resource allocation

