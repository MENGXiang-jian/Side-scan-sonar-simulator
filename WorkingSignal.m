%% WorkingSignal
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','FC','FR','BW','PW','DR','SP');
%% Original signal（TD）
fl = FC-0.5*BW;
fh = FC+0.5*BW;
fm = BW/PW;
pw = 0:1/FR:PW;
T = 0:1/FR:SP+15/c;

signal = real(exp(1j*2*pi*(fl*pw+fm*pw.*pw/2))); 
SIGNAL = zeros(1,length(T));
SIGNAL(1,1:length(signal)) = signal;
%% Original signal（FD）
SIGNAL_fft = 2*fft(SIGNAL);
TD_FFT = abs(SIGNAL_fft/length(SIGNAL));
FD = TD_FFT(1:length(SIGNAL)/2+1);
FD(2:end-1) = FD(2:end-1);
FD_F = FR*(0:(length(SIGNAL)/2))/length(SIGNAL);
%% Conversion signal(FD)
FL = 0; % fl-BW  0
FH = max(FD_F); % fh+BW  max(FD_F)
FLpoint = round((FL/FR)*length(SIGNAL))+1;
FHpoint = round((FH/FR)*length(SIGNAL))+1;
Narrowband = 0; % 0.01 0.0001
NarrowbandPoint = round(Narrowband*length(SIGNAL))+1;
SIGNAL_FPnum = FLpoint:NarrowbandPoint:FHpoint;   % Frequency Point Number
SIGNAL_FP = (SIGNAL_FPnum-1)*(FR)/length(SIGNAL); % Frequency Point
FP_F = SIGNAL_fft(SIGNAL_FPnum);                  % Frequency Point Fft
%% Conversion signal(TD)
FP_fft = zeros(1,length(SIGNAL_fft));
FP_fft(SIGNAL_FPnum) = FP_F;
SIGNAL_ifft = real(ifft(FP_fft));
signal_ifft = SIGNAL_ifft(1:length(pw));
%% Save data
save('WorkingSignal.mat','signal','SIGNAL','FR','T',...
    'signal_ifft','SIGNAL_ifft',...
    'SIGNAL_fft','SIGNAL_FP','SIGNAL_FPnum');
sound(sin(2*pi*10*(1:4000)/100));

%% Display
figure(1)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

% 1) Original signal（TD）
subplot(2,2,1);
plot(T,SIGNAL);
xlabel('s(t)');
set(gca,'FontSize',12);
xlim([0,length(pw)*(1/FR)]);
ylim([min(SIGNAL),max(SIGNAL)]);
% 2) Original signal（FD）
subplot(2,2,2);
plot(FD_F,FD);
xlabel('X(f)');
set(gca,'FontSize',12);
% 3) Conversion signal(TD)
subplot(2,2,3);
plot(T,SIGNAL_ifft);
xlabel('s"(t)');
set(gca,'FontSize',12);
xlim([0,length(pw)*(1/FR)]);
ylim([min(SIGNAL_ifft),max(SIGNAL_ifft)]);
% 4) Conversion signal(FD)
subplot(2,2,4);
plot(SIGNAL_FP,FD(SIGNAL_FPnum));
xlabel('X(f),f∈BW');
set(gca,'FontSize',12);

