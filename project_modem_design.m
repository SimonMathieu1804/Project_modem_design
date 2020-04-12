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

%need to vary Es/No to obtain BER in function of it
Es = 1; 
N0 = 1;

% vector of 2048 random bits
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

% Pulse Shapping + put at carrier freq ??

% AWGN Channel

% Matched filter and sampling + bring back from carrier freq ??

% 






