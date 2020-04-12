% LELEC2880: Modem design - Project
% Authors: DE COCK Justin, DELHAYE Quentin, SIMON Mathieu
% Date: 12/04/20

%%  Step 1 : Basic OFDM chain
%====== given values =========
N = 128; %number of subcarrier
f_0 = 2E9; %carrier frequency
f_sub = 15E3; %carrier subspacing
L = 16; %cyclic prefix length
%=============================
Nb = 2*N; % block size

%need to vary No to obtain BER in function of ES/N0
Es = 1; % keep Es to unity
N0 = 0;

% vector of 2048 random bits (to send 4 OFDM packets)
bits = randi([0 1],1,4*2*Nb); 

% Symbol mapping
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

% Seriel to parralel
parallel = [symbols(1:Nb) symbols(Nb+1:2*Nb) symbols(2*Nb+1:3*Nb) symbols(3*Nb+1:4*Nb)]; % each column is a block of 256 symbols

% IFFT on the blocks
parallel = ifft(parallel);

% Cyclic prefix insertion
CP = [parallel(end-L+1:end,1) parallel(end-L+1:end,2) parallel(end-L+1:end,3) parallel(end-L+1:end,4)];
paralel_CP = [CP ; parallel];

% parallel to serial
serial = [paralel_CP(:,1)' paralel_CP(:,2)' paralel_CP(:,3)' paralel_CP(:,4)'];

% Pulse Shapping 
alpha=0.2;
N_truncated=10;
u= rcosdesign(alpha,N_truncated,10,'sqrt'); 
E_u= u*u'; %E_u should be equal to one
fvtool(u,'impulse'); % Plot the filter
x = upfirdn(serial, u, 10);

% Shift at carrier freq ??

% AWGN Channel
x = x + randn(size(x))*2*N0;

% bring back from carrier freq ??

% Matched filter
y = upfirdn(x, u, 1, 10);

% serial to parralel 
y=y';
parallelRx = [y(1:(Nb+L)) y((Nb+L)+1:2*(Nb+L)) y(2*(Nb+L)+1:3*(Nb+L)) y(3*(Nb+L)+1:4*(Nb+L))];

% Remove CP
parallelRx = parallelRx(17:end,:);

% FFT on the blocks
parallelRx = fft(parallelRx);

% parallel to serial
output = [parallelRx(:,1)' parallelRx(:,2)' parallelRx(:,3)' parallelRx(:,4)'];
figure(3);
x = real(output); y = imag(output);
scatter(x,y,40,'o','filled','r'); title('Tx constellation','Fontsize',16); 
xlabel('In phase amplitude','Fontsize',14); ylabel('Quandrature amplitude','Fontsize',14);










