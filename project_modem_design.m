%% LELEC2880: Modem design - Project : Simulation
% Authors: DE COCK Justin, DELHAYE Quentin, SIMON Mathieu
% Date: 02/05/20

%% begining of the simulation
format long;
clear all;
clc;

%%  Step 1 : Basic OFDM chain
% 
% In this first step, the basic OFDM chain is implemented :
% 
% 1) Generate bits to send
% 2) Perform the symbol mapping (4-QAM)
% 3) Serial-to-Parallel computation (S/P)
% 4) Temporal representation of the blocks (IFFT)
% 5) Insert Cyclic Prefix (+CP)
% 6) Parallel-to-Serial (P/S)
% 7) Channel effect : AWGN (+real channel effect, Step3)
% 8) Serial-to-Parallel (S/P)
% 9) Remove Cyclic Prefix (-CP)
% 10) Frequency representation of the blocks (FFT)
% 11) Parallel-to-Serial (P/S)
% 12) Perform the symbol demapping 
% 13) Decision ($\hat{Bits}$)
% 
% Inputs :
%  - N (scalar) := Number of subcarrier
%  - L (scalar) := Cyclic Prefix length
%  - Nb (scalar) := Block Size
%  - N0 (scalar) := PSD of the noise
%  
% Outputs :
%   - Graphs showing the transmitted/received symbols 
%             and simulated/comparison BER curves
            
    


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
% Graph showing the transmitted constellation
figure(1);
x = real(symbols); y = imag(symbols);
scatter(x,y,40,'o','filled','r'); title('Tx constellation','Fontsize',16,'interpreter','latex');
xlabel('In phase amplitude','Fontsize',14,'interpreter','latex'); ylabel('Quadrature amplitude','Fontsize',14,'interpreter','latex');
% 3) Serial to parralel
parallel = [symbols(1:Nb) symbols(Nb+1:2*Nb) symbols(2*Nb+1:3*Nb) symbols(3*Nb+1:4*Nb)]; % each column is a block of 256 symbols
% 4) IFFT on the blocks
parallel = sqrt(Nb)*ifft(parallel);
% 5) Cyclic prefix insertion
CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4)];
paralel_CP = [CP ; parallel];
% 6) parallel to serial
serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).'];
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
% Graph showing the received constellation
figure(4);
x = real(output); y = imag(output);
scatter(x,y,40,'o','filled','r'); title('Rx constellation','Fontsize',16,'interpreter','latex');
xlabel('In phase amplitude','Fontsize',14,'interpreter','latex'); ylabel('Quadrature amplitude','Fontsize',14,'interpreter','latex');
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
        % 3) Serial to parralel
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
        % Simulated BER computation
        BER(index_SNR)=BER(index_SNR)+sum(output_bits.'~=bits);
    end
    N_symbols = 4*2*Nb*20;
    BER(index_SNR)=BER(index_SNR)/N_symbols;
end

%Theoritical BER computation
k=4;
M=2^k;
x=sqrt(3*k*Es_N0/(M-1));
theoretical_SER=erfc(sqrt(0.5*(10.^(Es_N0_dB/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(Es_N0_dB/10))))).^2;

% Graphs showing the datas
figure(5);
semilogy(Es_N0_dB,theoretical_SER/2,'-r','LineWidth',1.5);
hold on;
semilogy(Es_N0_dB,BER,'-xb','LineWidth',1.5,'MarkerSize',8);
grid;
xlabel('$E_S/N_0$ [dB]','interpreter','latex'); ylabel('BER','interpreter','latex'); legend('Theoretical (4QAM)','Simulated');
title('BER vs $E_S/N_0$','Fontsize',16,'interpreter','latex');

%% Step 2 : Resource allocation 
% In this step, the resource allocation is studied using two different
% criteria
%
% Algorithm of Water-Filling:
%1)Start with initial guess of mu (scalar, the "water-level")
%2)Compute corresponding powers and total required power
%3)Decrease or increase mu by some amount if required power is too large/small
%
%Also, another metric is used by minimizing the amount of errors (MSE minimization) :
%1)Start with initial guess of lambda (scalar)
%2)Compute corresponding powers and total required power
%3)Decrease or increase mu by some amount if required power is too large/small
%
%
% Inputs :
%  - N (scalar) := Number of subcarrier
%  - L (scalar) := Cyclic Prefix length
%  - Nb (scalar) := Block Size
%  - Pmax (scalar) := Maximum power available among all the subcarriers
%  
% Outputs :
%   - The powers allocated to each subcarrier ((Nb,1) vector)
%   - The bits/symbol, bits rate in each subcarrier ((Nb,1) vector)
%   - Graphs showing the channel frequency response, the powers allocated in each subcarrier,
%             the bits/symbol in each subcarrier.
            
N = 128; %number of subcarrier
L = 16; %cyclic prefix length
Nb = N; % block size

load('CIR.mat')
mu = 150; 
Pi = zeros(1,N);
hzeropad = [h.' zeros(1,N-length(h))];
Hf = fft(hzeropad);
Hf = fftshift(abs(Hf))/sqrt(8);

E_on_N=100;
Pmax=1;
sigmaN = 1/N/E_on_N;

figure(21);
plot(1:length(Hf),Hf,'LineWidth',1.5,'MarkerSize',8)
xlabel('frequency','interpreter','latex','Fontsize',14); ylabel('amplitude','interpreter','latex','Fontsize',14);
title('Channel frequency response','Fontsize',16,'interpreter','latex');

for (n=0:10000000)
    Pi = ((mu*(abs(Hf).^2)-sigmaN)./(abs(Hf).^2));
    Pi0 = Pi>0;
    Pi = Pi0.*Pi;
    Ptot = sum(Pi);
    if(abs(Pmax-Ptot)<0.1)
        %End of algorithm
        Ptot;
        break;
    elseif(Pmax>Ptot)
        %Too Few power
        mu = mu +0.001;
    elseif(Pmax<Ptot)
        %Too much power
        mu = mu - 0.001;
    end
end

f1 = figure(99);
clf;
set(f1,'Color',[1 1 1]);
bar(Pi +sigmaN./(abs(Hf).^2),1,'r')
hold on;
bar(sigmaN./abs(Hf).^2,1);
xlabel('subchannel indices','interpreter','latex','Fontsize',14);
title('Water filling algorithm','interpreter','latex','Fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%
%%%Bits performance%%%
%%%%%%%%%%%%%%%%%%%%%%
Petarg = 10^-5;
Gamma = 2/3*((erfcinv(Petarg/2))^2);
%- Water-filling distribution power
BitsWF = 1/2*log(1+((Pi.*(abs(Hf).^2))./(sigmaN*Gamma)))/log(2);%
%-Power uniformly distributed
BitsPowerUniform = 1/2*log(1+(Pmax/Nb.*(abs(Hf).^2).*ones(1,Nb)./(sigmaN*Gamma)))/log(2);

%Graph representation
figure();
hold on;
plot(1:128,BitsWF,'LineWidth',1.5,'MarkerSize',8);
plot(1:128,BitsPowerUniform,'LineWidth',1.5,'MarkerSize',8);
xlabel('Subcarrier','interpreter','latex','Fontsize',14); ylabel('bits per symbol','interpreter','latex','Fontsize',14); legend('Waterfilling','Uniform power allocation');
title('Number of bits per symbol for each subchannel','Fontsize',16,'interpreter','latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step2 :  Bonus : Bits fixed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = -0.5; %First Guess

for (n=0:1000000)
    Pi = sqrt(sigmaN./Hff.^2./(-lambda));
    Ptot = sum(Pi);
    if(abs(Pmax-Ptot)<0.01)
        %End of algorithm
        Ptot;
        break;
    elseif(Pmax>Ptot)
        %Too Few power
        lambda = lambda + 0.00001;
    elseif(Pmax<Ptot)
        %Too much power
        lambda = lambda - 0.00001;
    end
    %lambda
end

% Graph representation

f1 = figure();
clf;
set(f1,'Color',[1 1 1]);
bar(Pi +sigmaN./(abs(Hff).^2),1,'r')
hold on;
bar(sigmaN./abs(Hff).^2,1);
xlabel('subchannel indices','interpreter','latex','Fontsize',14);
title('Power allocation with MSE criterion','interpreter','latex','Fontsize',16)

Petarg = 10^-5;
Gamma = 2/3*((erfcinv(Petarg/2))^2);
%- MSE distribution power
BitsMSE = 1/2*log(1+(Pi.*(abs(Hf).^2))./(sigmaN*Gamma))/log(2);

figure();
hold on;
plot(1:128,BitsMSE,'LineWidth',1.5,'MarkerSize',8);
xlabel('Subcarrier','interpreter','latex','Fontsize',14); 
ylabel('bits per symbol','interpreter','latex','Fontsize',14); 
legend('MSE criterion allocation');
title('Number of bits per symbol for each subchannel','Fontsize',16,'interpreter','latex');



%% Step 3 : Channel estimation
% 
% In this step 3, the basic OFDM chain is extended by adding the channel estimation block :
% 1) Extraction of the received training sequence
% 2) Computation of the pseudo-inverse of the transmitted training sequence
% 3) Perform \hat{H} = Y*T^+
% 4) Temporal representation of the estimated channel
% 5) Using the prior knowledge, the extra-taps are zeroed out
% 
% Inputs :
%  - N (scalar) := Number of subcarrier
%  - L (scalar) := Cyclic Prefix length
%  - Nb (scalar) := Block Size
%  
% Outputs :
%   - Graphs showing the MSE between the real channel and the estimation
%                 in function of the SNR for different prior knowledges  

Nb = N;

% MSE calculation for each noise level and each length of channel
Nsnr=20; 
Es_N0_dB=linspace(0,20,Nsnr);
Es_N0=10.^(Es_N0_dB/10);
moy = 100;
for long = [4,7,8,9,12,120]
MSE=zeros(Nsnr,1);
for index_SNR=1:Nsnr
    
    for iter = 1:moy
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
        %bitstraining = [ones(1,Nb/2); -1*ones(1,Nb/2)]; %randi([0 1],1,Nb);
        maptraining = zeros(1,Nb);
        for i = 1:Nb/2
            maptraining(2*i-1) = 1;
            maptraining(2*i) = -1;
        end
        % 2) Symbol mapping
        training = [maptraining];
        % 3) Serial to parralel
        parallel = [training.' symbols(1:Nb) symbols(Nb+1:2*Nb) symbols(2*Nb+1:3*Nb) symbols(3*Nb+1:4*Nb)]; % each column is a block of 256 symbols
        % 4) IFFT on the blocks
        parallel = sqrt(Nb)*ifft(parallel);
        % 5) Cyclic prefix insertion
        CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4) parallel(end-L+1:end,5)];
        paralel_CP = [CP ; parallel];
        % 6) parallel to serial
        serial = [paralel_CP(:,1).' paralel_CP(:,2).' paralel_CP(:,3).' paralel_CP(:,4).' paralel_CP(:,5).'];
        % 7) AWGN channel + real channel
        %%%%%%%%%%%%%%%%%%%%%%
        %%%Channel creation%%%
        %%%%%%%%%%%%%%%%%%%%%%
        if (long==120)
            h = raylrnd(1:120) + 1i*raylrnd(1:120);
            h = h./norm(h);
            h = [h(1:120) zeros(1,(Nb-120))];
        else
            h = raylrnd(1:8) + 1i*raylrnd(1:8);
            h = h./norm(h);
            h = [h(1:8) zeros(1,(Nb-8))];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Channel impacts the data%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = conv(h,serial)+ randn(size(conv(h,serial)))*sqrt(N0/2)+ randn(size(conv(h,serial)))*sqrt(N0/2)*1i;
        % 8) serial to parralel
        y=y.';
        parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L)) y(4*(Nb+L)+1:5*(Nb+L))];
        % 9) Remove CP
        parallelRx = parallelRx((L+1):end,:);
        % 10) FFT on the blocks
        parallelRx = fft(parallelRx)/sqrt(Nb);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Channel estimation%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Extraction of the received training sequence
        outtraining = parallelRx(:,1);
        % Computation of the pseudo-inverse of the 
        %    transmitted training sequence
        traindague = pinv(training, 10^-5);
        % Perform \hat{H} = Y*T^+
        hhatfreq = traindague.*outtraining;
        % Temporal representation of the estimated channel
        hguess = (Nb)*(ifft(hhatfreq));
        % Using the prior knowledge, the extra-taps are zeroed out
        estimee = [hguess(1:long).' zeros(1,(Nb-long))];
        
        %%%%%%%%%%%%%%%%%%%%%
        %%%MSE computation%%%
        %%%%%%%%%%%%%%%%%%%%%
        MSECalc = sum(abs(estimee-h).^2)/moy;
        MSE(index_SNR) = MSE(index_SNR) + MSECalc;
        
        %%%%%%%%%%%%%%%%%%%
        %%% End of Step3%%%
        %%%%%%%%%%%%%%%%%%%
        
        % 11) parallel to serial
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
    end
end
%Graph showing the MSE in function of the SNR
figure(100);
plot(Es_N0_dB,10*log10(MSE),'-x','linewidth',1.5);
title('The mean-square error of the channel estimation in function of the SNR','interpreter','latex','Fontsize',16);
xlabel('SNR [dB]','interpreter','Latex','Fontsize',14);
ylabel('MSE [dB]','interpreter','latex','Fontsize',14);
hold on;
legend('4 taps','7 taps','8 taps','9 taps','12 taps','120 taps');
end

%% Step 4 : Optimal Viterbi (soft) decoding
% 
% In this step 4, the basic OFDM chain is extended by adding the
% equalization and the Viterbi encoding/decoding blocks.
% Two different cases are studied :
% 1) Viterbi decoding with an AWGN channel
% 2) Viterbi decoding with a real +AWGN channel using the real or estimated
% channel for the equalization or Viterbi decoding
% 
% Inputs :
%  - N (scalar) := Number of subcarrier
%  - L (scalar) := Cyclic Prefix length
%  - Nb (scalar) := Block Size
%  - g (matrix) := Convolutional code
%
% Outputs :
%   - Graphs showing the BER 
%                 in function of the SNR for different types of channel  

%====================== With an AWGN channel ==============================

N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
Nb = N; % block size
L = 16; % Cyclic prefix length
Lf = 128; % Input sequence size
%Lf = 8; % Test input sequence size
g = [1 1 0;1 1 1]; % Convolutional code
[n,K] = size(g); % length of the code
k = 1; % rank of the code
N0 = 0.01;

% 1) generate 128 random bits and code them
u = randi([0 1],1,Lf); % Input sequence
x = zeros(n,Lf); % Coded sequence
u_pad = [zeros(1,n) u];
for i = 1:Lf
    t = i+2;
    x(1,i) = sum(u_pad(i:t).*(flip(g(1,:))));
    x(2,i) = sum(u_pad(i:t).*(flip(g(2,:))));
end
x(x==3) = 1;
x(x==2) = 0;
% 2) Symbol mapping
map = x;
map(map==0) = -1;
map = sqrt(2)/2*map;
symbols = zeros(Lf,1); % One OFDM symbol
for r=1:Lf
    symbols(r)=map(1,r)+1i*map(2,r);
end
% 3) Seriel to parralel
parallel = symbols; % each column is a block of 128 symbols
% 4) IFFT on the blocks
parallel = sqrt(Nb)*ifft(parallel);
% 5) Cyclic prefix insertion
CP = parallel(end-L+1:end);
paralel_CP = [CP ; parallel];
% 6) parallel to serial
serial = paralel_CP(:,1).';
% 7) AWGN channel
y = serial+ randn(size(serial))*sqrt(N0/2)+ randn(size(serial))*sqrt(N0/2)*1i;
% 8) serial to parralel
y=y.';
parallelRx = y(1:(Lf+L));
% 9) Remove CP
parallelRx = parallelRx((L+1):end,:);
% 10) FFT on the blocks
parallelRx = fft(parallelRx)/sqrt(Nb);
% 11) parallel to serial
output = parallelRx(:,1).';
% 12) demapping
coded_output_bits = zeros(Lf,1);
for r=1:Lf
    coded_output_bits(2*r-1)=real(output(r));
    coded_output_bits(2*r)=imag(output(r));
end

%---------------------------- Viterbi decoding

% states are 0, 1, 2 or 3
state_table=[0]; % starting state
distance = [0]; % starting distance

% find the best path
for iter = 1:Lf
    
    % create the matrix of the new possible paths
    a1 = state_table(:,end);
    a1(a1==0)= 0;
    a1(a1==2)= 0;
    a1(a1==1)= 2;
    a1(a1==3)= 2;
    a2 = state_table(:,end);
    a2(a2==3)= 3;
    a2(a2==1)= 3;
    a2(a2==0)= 1;
    a2(a2==2)= 1;
    a = [a1 ;a2];
    state_table = [state_table ; state_table];
    state_table = [state_table a];
    
    % update the matrix of cumulative distance
    distance = [distance distance];
    twolast = [state_table(:,end-1) state_table(:,end)];
    [not, yep] = size(distance);
    for k = 1:yep
        switch(twolast(k,1))
            case 0
                if(twolast(k,2)==0)
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                else
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                end
            case 1
                if(twolast(k,2)==2)
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                else
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                end
            case 2
                if(twolast(k,2)==0)
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                else
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                end
            case 3
                if(twolast(k,2)==2)
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                else
                    distance(1,k)=distance(1,k)+(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                end
            otherwise
                fprintf('error');
        end
    end
    
    % supress paths that lead to a same state but with higher cumuative distance
    result = [];
    newDist = [];
    for k = 0:3
        index = find(state_table(:,end)==k);
        if (isempty(index)==0)
            A = state_table(index,:);
            B = distance(1,index);
            [M I] = min(B);
            result = [result ; A(I,:)];
            newDist = [newDist B(I)];
        end
    end
    state_table = result;
    distance = newDist;
end
[best_dist, indpath] = min(distance);
best_path = state_table(indpath,:);

% interpreting the best path into the bit sequence
output = mod(best_path(2:end),2);
error = sum(abs(u - output)); 


%====================== With a channel h ==================================

N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
Nb = N; % block size
L = 16; % Cyclic prefix length
Lf = 128; % Input sequence size
%Lf = 8; % Test input sequence size
g = [1 1 0;1 1 1]; % Convolutional code
[n,K] = size(g); % length of the code
k = 1; % rank of the code

Nber=20; %%Changed 20 before
Es_N0_dB=linspace(0,20,Nber);
Es_N0=10.^(Es_N0_dB/10);
BER1=zeros(Nber,1); % to make the graph for Modified Viterbi with estimation
BER2=zeros(Nber,1); % to make the graph for Classical Viterbi
BER3=zeros(Nber,1); % to make the graph for Modified Viterbi with perfect knowledge of h
BER4=zeros(Nber,1); %to make the graph for no coding

% run several time in order to make the different graphs
for graph = 1:4
    
    BER=zeros(Nber,1);
    for index_BER=1:Nber
        for iterations = 1:50
            N0=1/Es_N0(index_BER);
            
            % 1) vector of 128 random bits and code it
            u = randi([0 1],1,Lf); 
            x = zeros(n,Lf); % Coded sequence
            u_pad = [zeros(1,n) u];
            for i = 1:Lf
                t = i+2;
                x(1,i) = sum(u_pad(i:t).*(flip(g(1,:))));
                x(2,i) = sum(u_pad(i:t).*(flip(g(2,:))));
            end
            x(x==3) = 1;
            x(x==2) = 0;
            if(graph==4)
                x(1,:)= randi([0 1],1,Lf);
                x(2,:)= randi([0 1],1,Lf);
            end
            % 2) Symbol mapping and insert the training sequence
            map = x;
            map(map==0) = -1;
            map = sqrt(2)/2*map;
            symbols = zeros(Lf,1); % One OFDM symbol
            for r=1:Lf
                symbols(r)=map(1,r)+1i*map(2,r);
            end
            bitstraining = [ones(1,Lf/2); -1*ones(1,Lf/2)]; %randi([0 1],1,Nb);
            maptraining = zeros(1,Lf);
            for i = 1:Lf/2
                maptraining(2*i-1) = bitstraining(1,i);
                maptraining(2*i) = bitstraining(2,i);
            end
            training = [maptraining];
            % 3) Seriel to parralel
            parallel = [training.' symbols(1:Nb)];
            % 4) IFFT on the blocks
            parallel = sqrt(Lf)*ifft(parallel);
            % 5) Cyclic prefix insertion
            CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2)];
            paralel_CP = [CP ; parallel];
            % 6) parallel to serial
            serial = [paralel_CP(:,1).' paralel_CP(:,2).'];
            % 7) Create the real channel and apply it on the data
            long = 8;
            h = raylrnd(1:long) + 1i*raylrnd(1:long);
            h = h./norm(h);
            y = conv(h,serial)+ randn(size(conv(h,serial)))*sqrt(N0/2)+ randn(size(conv(h,serial)))*sqrt(N0/2)*1i;
            % 8) serial to parralel
            y=y.';
            parallelRx = [y(1:(Lf+L)) y((Lf+L)+1:2*(Lf+L))];
            % 9) Remove CP
            parallelRx = parallelRx((L+1):end,:);
            % 10) FFT on the blocks
            parallelRx = fft(parallelRx)/sqrt(Lf);
            % 11) channel estimation
            test = parallelRx(:,1);
            testtraining = test;
            traindague = pinv(training, 10^-5);
            hhat = traindague.*testtraining;
            hguess = (Lf)*(ifft(hhat));
            estimee = [hguess(1) hguess(2) hguess(3) hguess(4) hguess(5) hguess(6) hguess(7) hguess(8) zeros(1,Lf-8)];
            if(graph==3)
                estimee = [h zeros(1,Lf-8)];
            end
            estimation = fft(estimee)/Lf;
            % 12) parallel to serial and equalization
            outtraining = parallelRx(:,1).';
            output_seq = parallelRx(:,2).'./estimation;
            % 13) demapping
            output_bits = zeros(Lf,1);
            for k=1:Lf
                output_bits(2*k-1)=real(output_seq(k));
                output_bits(2*k)=imag(output_seq(k));
            end
            coded_output_bits = output_bits;
            if(graph==2)
                estimation = ones(1,Lf);
            end
            
            %--------------------------- Viterbi decoding 
            
            % states are 0, 1, 2 or 3
            state_table=[0]; % starting state
            distance = [0]; % starting distance
            
            % find the best path
            for iter = 1:Lf
                
                % create the matrix of the new possible paths
                a1 = state_table(:,end);
                a1(a1==0)= 0;
                a1(a1==2)= 0;
                a1(a1==1)= 2;
                a1(a1==3)= 2;
                a2 = state_table(:,end);
                a2(a2==3)= 3;
                a2(a2==1)= 3;
                a2(a2==0)= 1;
                a2(a2==2)= 1;
                a = [a1 ;a2];
                state_table = [state_table ; state_table];
                state_table = [state_table a];
                
                % update the matrix of cumulative distance
                distance = [distance distance];
                twolast = [state_table(:,end-1) state_table(:,end)];
                [not, yep] = size(distance);
                for k = 1:yep
                    switch(twolast(k,1))
                        case 0
                            if(twolast(k,2)==0)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            end
                        case 1
                            if(twolast(k,2)==2)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            end
                        case 2
                            if(twolast(k,2)==0)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            end
                        case 3
                            if(twolast(k,2)==2)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            end
                        otherwise
                            fprintf('error');
                    end
                end
                
                % supress paths that lead to a same state but with higher cumuative distance
                result = [];
                newDist = [];
                for k = 0:3
                    index = find(state_table(:,end)==k);
                    if (isempty(index)==0)
                        A = state_table(index,:);
                        B = distance(1,index);
                        [M I] = min(B);
                        result = [result ; A(I,:)];
                        newDist = [newDist B(I)];
                    end
                end
                state_table = result;
                distance = newDist;
            end
            
            [best_dist, indpath] = min(distance);
            best_path = state_table(indpath,:);
            
            % interpreting the best path into the bit sequence
            output = mod(best_path(2:end),2);
            
            %---------------------------------
            if(graph==4)
                output = coded_output_bits;
                output(output>0)=1;
                output(output<0)=0;
                fourth = [x(1,:); x(2,:)];
                fourth = fourth(:).';
                fourth = fourth.';
                BER(index_BER)=BER(index_BER)+sum(output~=fourth);
            else
                BER(index_BER)=BER(index_BER)+sum(output~=u);
            end
        end
        N_bits = Lf*50;
        if(graph==4)
            N_bits=N_bits*2;
        end
        BER(index_BER)=BER(index_BER)/N_bits;
    end
    if(graph == 1)
        BER1=BER;
    elseif (graph == 2)
        BER2=BER;
    elseif (graph ==3)
        BER3=BER;
    else 
        BER4=BER;
    end
end

semilogy(Es_N0_dB,BER2,'-xb','LineWidth',1.5,'MarkerSize',8);
hold on; 
semilogy(Es_N0_dB,BER1,'-xr','LineWidth',1.5,'MarkerSize',8);
hold on;
semilogy(Es_N0_dB,BER3,'-xg','LineWidth',1.5,'MarkerSize',8);
hold on;
semilogy(Es_N0_dB,BER4,'-xm','LineWidth',1.5,'MarkerSize',8);
grid;
xlabel('$E_S/N_0$ [dB]','Fontsize',14,'interpreter','latex'); ylabel('BER','Fontsize',14,'interpreter','latex'); 
legend('Classical Viterbi','Modified Viterbi with estimation','Modified Viterbi with perfect knowledge', 'without coding');
title('Viterbi: BER vs $E_S/N_0$','Fontsize',16,'interpreter','latex');



%% Step 4 Bonus
% 
% In this step 4, the basic OFDM chain is extended by adding the
% equalization and the Viterbi encoding/decoding blocks.
% The interleaving impact of the bits is studied.
% 
% Inputs :
%  - N (scalar) := Number of subcarrier
%  - L (scalar) := Cyclic Prefix length
%  - Nb (scalar) := Block Size
%  - g (matrix) := Convolutional code
% Outputs :
%   - Graphs showing the BER
%                 in function of the SNR with/without interleaving  

N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
Nb = N; % block size
L = 16; % Cyclic prefix length
Lf = 128; % Input sequence size
%Lf = 8; % Test input sequence size
g = [1 1 0;1 1 1]; % Convolutional code
[n,K] = size(g); % length of the code
k = 1; % rank of the code

Nber=20; %%Changed 20 before
Es_N0_dB=linspace(0,20,Nber);
Es_N0=10.^(Es_N0_dB/10);
BER1=zeros(Nber,1); % to make the graph for Modified Viterbi with estimation
BER2=zeros(Nber,1); % to make the graph for Classical Viterbi
BER3=zeros(Nber,1); % to make the graph for Modified Viterbi with perfect knowledge of h
BER4=zeros(Nber,1); %to make the graph for no coding

long = 8;
iter1 = 1000;
h_matrix = zeros(Nber,long,iter1);

random_matrix = linspace(1,Lf,Lf);
random_matrix = random_matrix(randperm(length(random_matrix)));

for index_BER = 1:Nber
    for iterations = 1:iter1
        h = raylrnd(1:long) + 1i*raylrnd(1:long);
        h = h./norm(h);
        h_matrix(index_BER,:,iterations) = h;
    end
end
% run several time in order to make the different graphs
for graph = 1:5
    
    BER=zeros(Nber,1);
    for index_BER=1:Nber
        for iterations = 1:iter1
            N0=1/Es_N0(index_BER);
            
            % 1) vector of 128 random bits and code it
            u = randi([0 1],1,Lf); 
            x = zeros(n,Lf); % Coded sequence
            u_pad = [zeros(1,n) u];
            for i = 1:Lf
                t = i+2;
                x(1,i) = sum(u_pad(i:t).*(flip(g(1,:))));
                x(2,i) = sum(u_pad(i:t).*(flip(g(2,:))));
            end
            x(x==3) = 1;
            x(x==2) = 0;
            x2_rand = zeros(1,Lf);
            for r1=1:Lf
                rand_ind1 = random_matrix(r1);
                x2_rand(r1) = x(2,rand_ind1);
            end
            if(graph==5)
                x(2,:) = x2_rand;
            end
            if(graph==4)
                x(1,:)= randi([0 1],1,Lf);
                x(2,:)= randi([0 1],1,Lf);
            end
            % 2) Symbol mapping and insert the training sequence
            map = x;
            map(map==0) = -1;
            map = sqrt(2)/2*map;
            symbols = zeros(Lf,1); % One OFDM symbol
            for r=1:Lf
                symbols(r)=map(1,r)+1i*map(2,r);
            end
            bitstraining = [ones(1,Lf/2); -1*ones(1,Lf/2)]; %randi([0 1],1,Nb);
            maptraining = zeros(1,Lf);
            for i = 1:Lf/2
                maptraining(2*i-1) = bitstraining(1,i);
                maptraining(2*i) = bitstraining(2,i);
            end
            training = [maptraining];
            % 3) Seriel to parralel
            parallel = [training.' symbols(1:Nb)];
            % 4) IFFT on the blocks
            parallel = sqrt(Lf)*ifft(parallel);
            % 5) Cyclic prefix insertion
            CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2)];
            paralel_CP = [CP ; parallel];
            % 6) parallel to serial
            serial = [paralel_CP(:,1).' paralel_CP(:,2).'];
            % 7) Create the real channel and apply it on the data
%             long = 8;
%             h = raylrnd(1:long) + 1i*raylrnd(1:long);
%             h = h./norm(h);
            h = h_matrix(index_BER,:,iterations);
            y = conv(h,serial)+ randn(size(conv(h,serial)))*sqrt(N0/2)+ randn(size(conv(h,serial)))*sqrt(N0/2)*1i;
            % 8) serial to parralel
            y=y.';
            parallelRx = [y(1:(Lf+L)) y((Lf+L)+1:2*(Lf+L))];
            % 9) Remove CP
            parallelRx = parallelRx((L+1):end,:);
            % 10) FFT on the blocks
            parallelRx = fft(parallelRx)/sqrt(Lf);
            % 11) channel estimation
            test = parallelRx(:,1);
            testtraining = test;
            traindague = pinv(training, 10^-5);
            hhat = traindague.*testtraining;
            hguess = (Lf)*(ifft(hhat));
            estimee = [hguess(1) hguess(2) hguess(3) hguess(4) hguess(5) hguess(6) hguess(7) hguess(8) zeros(1,Lf-8)];
            if(graph==3 || graph==5)
                estimee = [h zeros(1,Lf-8)];
            end
            estimation = fft(estimee)/Lf;
            % 12) parallel to serial and equalization
            outtraining = parallelRx(:,1).';
            output_seq = parallelRx(:,2).'./estimation;
            % 13) demapping
            output_bits = zeros(Lf,1);
            x2_not_rand = zeros(Lf,1);
            for r2=1:Lf
                rand_ind2 = random_matrix(r2);
                x2_not_rand(rand_ind2) = imag(output_seq(r2));
            end
            for k=1:Lf
                output_bits(2*k-1)=real(output_seq(k));
                if (graph == 5)
                    output_bits(2*k) = x2_not_rand(k);
                else
                    output_bits(2*k)=imag(output_seq(k));
                end
            end
            coded_output_bits = output_bits;
            if(graph==2)
                estimation = ones(1,Lf);
            end
            
            %--------------------------- Viterbi decoding 
            
            % states are 0, 1, 2 or 3
            state_table=[0]; % starting state
            distance = [0]; % starting distance
            
            % find the best path
            for iter = 1:Lf
                
                % create the matrix of the new possible paths
                a1 = state_table(:,end);
                a1(a1==0)= 0;
                a1(a1==2)= 0;
                a1(a1==1)= 2;
                a1(a1==3)= 2;
                a2 = state_table(:,end);
                a2(a2==3)= 3;
                a2(a2==1)= 3;
                a2(a2==0)= 1;
                a2(a2==2)= 1;
                a = [a1 ;a2];
                state_table = [state_table ; state_table];
                state_table = [state_table a];
                
                % update the matrix of cumulative distance
                distance = [distance distance];
                twolast = [state_table(:,end-1) state_table(:,end)];
                [not, yep] = size(distance);
                for k = 1:yep
                    if(graph==5)
                    switch(twolast(k,1))
                        case 0
                            if(twolast(k,2)==0)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                            end
                        case 1
                            if(twolast(k,2)==2)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                            end
                        case 2
                            if(twolast(k,2)==0)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                            end
                        case 3
                            if(twolast(k,2)==2)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)-sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)+sqrt(2)/2)^2;
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*(coded_output_bits(2*iter-1)+sqrt(2)/2)^2+abs(estimation(find(random_matrix==iter)))^2*(coded_output_bits(2*iter)-sqrt(2)/2)^2;
                            end
                        otherwise
                            fprintf('error');
                    end
                    else
                        switch(twolast(k,1))
                        case 0
                            if(twolast(k,2)==0)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            end
                        case 1
                            if(twolast(k,2)==2)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            end
                        case 2
                            if(twolast(k,2)==0)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            end
                        case 3
                            if(twolast(k,2)==2)
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)-sqrt(2)/2)^2+(coded_output_bits(2*iter)+sqrt(2)/2)^2);
                            else
                                distance(1,k)=distance(1,k)+abs(estimation(iter))^2*((coded_output_bits(2*iter-1)+sqrt(2)/2)^2+(coded_output_bits(2*iter)-sqrt(2)/2)^2);
                            end
                        otherwise
                            fprintf('error');
                    end
                    end
                end
                
                % supress paths that lead to a same state but with higher cumuative distance
                result = [];
                newDist = [];
                for k = 0:3
                    index = find(state_table(:,end)==k);
                    if (isempty(index)==0)
                        A = state_table(index,:);
                        B = distance(1,index);
                        [M I] = min(B);
                        result = [result ; A(I,:)];
                        newDist = [newDist B(I)];
                    end
                end
                state_table = result;
                distance = newDist;
            end
            
            [best_dist, indpath] = min(distance);
            best_path = state_table(indpath,:);
            
            % interpreting the best path into the bit sequence
            output = mod(best_path(2:end),2);
            
            %---------------------------------
            if(graph==4)
                output = coded_output_bits;
                output(output>0)=1;
                output(output<0)=0;
                fourth = [x(1,:); x(2,:)];
                fourth = fourth(:).';
                fourth = fourth.';
                BER(index_BER)=BER(index_BER)+sum(output~=fourth);
            else
                BER(index_BER)=BER(index_BER)+sum(output~=u);
            end
        end
        N_bits = Lf*iter1;
        if(graph==4)
            N_bits=N_bits*2;
        end
        BER(index_BER)=BER(index_BER)/N_bits;
    end
    if(graph == 1)
        BER1=BER;
    elseif (graph == 2)
        BER2=BER;
    elseif (graph ==3)
        BER3=BER;
    elseif (graph==4)
        BER4=BER;
    elseif (graph==5)
        BER5 = BER;
    end
end

% semilogy(Es_N0_dB,BER2,'-xb','LineWidth',1.5,'MarkerSize',8);
% hold on; 
% semilogy(Es_N0_dB,BER1,'-xr','LineWidth',1.5,'MarkerSize',8);
% hold on;
figure();
semilogy(Es_N0_dB,BER3,'-xg','LineWidth',1.5,'MarkerSize',8);
hold on;
% semilogy(Es_N0_dB,BER4,'-xm','LineWidth',1.5,'MarkerSize',8);
% hold on;
semilogy(Es_N0_dB,BER5,'-x','LineWidth',1.5,'MarkerSize',8);
grid;
xlabel('$E_S/N_0$ [dB]','Fontsize',14,'interpreter','latex'); ylabel('BER','Fontsize',14,'interpreter','latex'); 
% legend('Classical Viterbi','Modified Viterbi with estimation','Modified Viterbi with perfect knowledge', 'without coding');
legend('Without interleaving','With interleaving');
title('Viterbi: BER vs $E_S/N_0$','Fontsize',16,'interpreter','latex');
