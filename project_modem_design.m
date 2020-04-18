%% LELEC2880: Modem design - Project : Simulation
% Authors: DE COCK Justin, DELHAYE Quentin, SIMON Mathieu
% Date: 12/04/20

%% begining of the simulation
format long;
clear all;
clc;

%%  Step 1 : Basic OFDM chain
%======================= given values =====================================
N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
Nb = N; % block size

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
parallel = sqrt(Nb)*ifft(parallel);
% 5) Cyclic prefix insertion
CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4)];
paralel_CP = [CP ; parallel];
% 6) parallel to serial
serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).'];

% % Pulse Shapping
% alpha=0.2;
% N_truncated=10;
% u= rcosdesign(alpha,N_truncated,10,'sqrt');
% E_u= u*u'; %E_u should be equal to one
% fvtool(u,'impulse'); % Plot the filter
% fvtool(u,'freq');
% x = upfirdn(serial, u, 10);
% % Shift at carrier freq ??
% % AWGN Channel
% x = x + randn(size(x))*2*N0;
% % bring back from carrier freq ??
% % Matched filter
% y = upfirdn(x, u, 1, 10);

% 7) AWGN channel
y = serial+ randn(size(serial))*sqrt(N0/2)+ randn(size(serial))*sqrt(N0/2)*1i; % line to comment in order to use pulse shaping
% 8) serial to parralel
y=y.';
parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L))];
% 9) Remove CP
parallelRx = parallelRx((L+1):end,:);
% 10) FFT on the blocks
parallelRx = fft(parallelRx)/sqrt(Nb);
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

BER = sum(output_bits.'~=bits)/2048;

%======================= Es=1 and varying N0 ==============================

% BER calculation for each noise level
Nsnr=20;
Es_N0_dB=linspace(0,12,Nsnr);
Es_N0=10.^(Es_N0_dB/10);

BER=zeros(Nsnr,1);
for index_SNR=1:Nsnr
    for iter = 1:20
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
        parallel = sqrt(Nb)*ifft(parallel);
        % 5) Cyclic prefix insertion
        CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4)];
        paralel_CP = [CP ; parallel];
        % 6) parallel to serial
        serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).'];
        % 7) AWGN channel
        y = serial+ randn(size(serial))*sqrt(N0/2)+ randn(size(serial))*sqrt(N0/2)*1i;
        % 8) serial to parralel
        y=y.';
        parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L))];
        % 9) Remove CP
        parallelRx = parallelRx((L+1):end,:);
        % 10) FFT on the blocks
        parallelRx = fft(parallelRx)/sqrt(Nb);
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
        %N_symbols = 4*2*Nb;
        BER(index_SNR)=BER(index_SNR)+sum(output_bits.'~=bits);%/N_symbols;
    end
    N_symbols = 4*2*Nb*5;
    BER(index_SNR)=BER(index_SNR)/N_symbols;
end


k=4;
M=2^k;
x=sqrt(3*k*Es_N0/(M-1));
theoretical_BER=erfc(sqrt(0.5*(10.^(Es_N0_dB/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(Es_N0_dB/10))))).^2;

figure(5);
semilogy(Es_N0_dB,theoretical_BER,'-r','LineWidth',1.5);
hold on;
semilogy(Es_N0_dB,BER/2,'-xb','LineWidth',1.5,'MarkerSize',8);
grid;
xlabel('E_S/N_0 [dB]'); ylabel('SER'); legend('Theory (4QAM)','Simulated');

%% Step 2 bonus : Power allocation only

%First, computation Water-Filling
%1)Start with initial guess of �  =1/(2*lambda*ln(2))
%2)Compute corresponding powers and total required power
%3)Decrease or increase � by some amount if required power is too large/small
N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
Nb = N; % block size

load('CIR.mat')
%fvtool(h,'impulse'); % Plot the filter
%fvtool(h,'freq');
Pmax = 128; %1
Ptot = 0;
mu = 0.5; %0.5
Pi = zeros(1,N);
N0 = 0.01;%0.000035;
size(h)
size(zeros(1,N-length(h)))
%h = h/norm(h);
hzeropad = h;
hzeropad = [h.' zeros(1,N-length(h))];
Hf = fft(hzeropad);
Hf = fftshift(abs(Hf));%/sqrt(8);
%Hf = Hf/norm(Hf);
figure(21);
plot(1:length(Hf),Hf);

%Hf = [h(1)*ones(1,16) h(2)*ones(1,16) h(3)*ones(1,16) h(4)*ones(1,16) h(5)*ones(1,16) h(6)*ones(1,16) h(7)*ones(1,16) h(8)*ones(1,16)];

for (n=0:1000000)
    Pi = (mu*(abs(Hf).^2)-N0)./(abs(Hf).^2);
    Pi0 = Pi>0;
    Pi = Pi0.*Pi;
    Ptot = sum(Pi);
    if (abs(Pmax-Ptot)<0.01)
     %End of algorithm
        Ptot;
        break; 
    elseif(Pmax > Ptot)
        %Too Few power 
        mu = mu +0.0001;
    elseif (Pmax<Ptot)
        %Too much power
        mu = mu - 0.0001;
    end 
end

%Hff = [h(1)*ones(1,16) h(2)*ones(1,16) h(3)*ones(1,16) h(4)*ones(1,16) h(5)*ones(1,16) h(6)*ones(1,16) h(7)*ones(1,16) h(8)*ones(1,16)];
Hff = Hf;

f1 = figure();
    clf;
    set(f1,'Color',[1 1 1]);
    bar(Pi +N0./(abs(Hff).^2),1,'r')
    hold on;    
    bar(N0./abs(Hff).^2,1);
    xlabel('subchannel indices');
    title('Water filling algorithm')
    
   % legend('amount of power allocated to each subchannel',...
         %  'Noise to Carrier Ratio')
%%%%%%%%%%%%%%%%%%%%%%
%%%Bits performance
%%%%%%%%%%%%%%%%%%%%%%
Hff = Hff;
Petarg = 10^-5;
Gamma = 2/3*((erfcinv(Petarg/2))^2);
%- Water-filling distribution power
BitsWF = 1/2*log(1+(Pi.*(abs(Hff).^2))./(N0*Gamma))/log(2);%
%-Power uniformly distributed
BitsPowerUniform = 1/2*log(1+(Pmax/Nb*ones(1,Nb).*(abs(Hff).^2))./(N0*Gamma))/log(2);%.*(abs(Hff).^2 

figure();
hold on;
plot(1:128,BitsWF);
plot(1:128,BitsPowerUniform);
legend('WF','Uniform');

%%%%%%%%%%%%%%%%%%%%
% Step2 : False Bonus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Step 3 : Channel estimation
Nb = N;
% BER calculation for each noise level
Nsnr=1; %%Changed 20 before
Es_N0_dB=linspace(0,12,Nsnr);
Es_N0=10.^(Es_N0_dB/10);

BER=zeros(Nsnr,1);
for index_SNR=1:Nsnr
    for iter = 1:20
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
        bitstraining = randi([0 1],1,Nb);
        % 2) Symbol mapping
        maptraining = bitstraining;
        maptraining(maptraining==0) = -1;
        %map = sqrt(2)/2*map;
        training = [maptraining];
        % 3) Seriel to parralel
        size(training)
        size(symbols(1:Nb))
        parallel = [training.' symbols(1:Nb) symbols(Nb+1:2*Nb) symbols(2*Nb+1:3*Nb) symbols(3*Nb+1:4*Nb)]; % each column is a block of 256 symbols
        % 4) IFFT on the blocks
        parallel = sqrt(Nb)*ifft(parallel);
        % 5) Cyclic prefix insertion
        CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4) parallel(end-L+1:end,5)];
        paralel_CP = [CP ; parallel];
        % 6) parallel to serial
        serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).' paralel_CP(:,5).'];
        %y = serial+ randn(size(serial))*sqrt(N0/2)+ randn(size(serial))*sqrt(N0/2)*1i;
        % 7) AWGN channel + real channel
        %h = +0*1i;
        long = 8;
        h = raylrnd(1:long) + 1i*raylrnd(1:long);
        for i= 1:999
            h = raylrnd(1:long) + 1i*raylrnd(1:long);
        end
        h = h/1000;
        %h = raylrnd(1:8) + 1i*raylrnd(1:8);
        h = h./norm(h);
        y = conv(h,serial)+ randn(size(conv(h,serial)))*sqrt(N0/2)+ randn(size(conv(h,serial)))*sqrt(N0/2)*1i;
        size(h)
        size(serial)
        % 8) serial to parralel
        y=y.';
        parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L)) y(4*(Nb+L)+1:5*(Nb+L))];
        % 9) Remove CP
        parallelRx = parallelRx((L+1):end,:);
        % 10) FFT on the blocks
        parallelRx = fft(parallelRx)/sqrt(Nb);
        % CHANNEL ESTIMATION
        test = parallelRx(:,1);
        %testtraining = test((L+1):end);
        testsize = size(test)
        testtraining = test;
        traindague = pinv(training, 10^-5);
        %training = training.',
        size(traindague)
        size(training)
       % identity = traindague.'.*training;
        size(traindague)
        size(testtraining)
        hhat = traindague.*testtraining;
        %hhat = hhat.';
      %  figure(10);
       % plot(1:8,[mean(abs(hhat(1,:))) mean(abs(hhat(2,:))) mean(abs(hhat(3,:))) mean(abs(hhat(4,:))) mean(abs(hhat(5,:))) mean(abs(hhat(6,:))) mean(abs(hhat(7,:))) mean(abs(hhat(8,:)))]);
        hhat = hhat(1:16:end) % taking Each 16 elements of Hf in order to have 8 taps (Sampling method)
        figure(11);
        plot(1:length(h),abs(h));
        %estimee = [mean(abs(hhat(1,:))) mean(abs(hhat(2,:))) mean(abs(hhat(3,:))) mean(abs(hhat(4,:))) mean(abs(hhat(5,:))) mean(abs(hhat(6,:))) mean(abs(hhat(7,:))) mean(abs(hhat(8,:)))]
        hguess = (ifft(hhat))
        hguess = hguess/norm(hguess);
        estimee = [hguess(1) hguess(2) hguess(3) hguess(4) hguess(5) hguess(6) hguess(7) hguess(8)]
        (h)
        figure(20);
        plot(1:length(hguess), abs(hguess))
        size(hhat)
        %abs(hhat(2,:))
        % 11) parallel to serial
        outtraining = parallelRx(:,1).';
        output = [parallelRx(:,2).' parallelRx(:,3).' parallelRx(:,4).' parallelRx(:,5).'];
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
        %N_symbols = 4*2*Nb;
        BER(index_SNR)=BER(index_SNR)+sum(output_bits.'~=bits);%/N_symbols;
    end
    N_symbols = 4*2*Nb*5;
    BER(index_SNR)=BER(index_SNR)/N_symbols;
end
