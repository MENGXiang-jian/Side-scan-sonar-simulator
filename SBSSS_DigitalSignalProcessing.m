%% SBSSS Digital signal processing
%% Tips


%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','TRANSMITTER','FC','BW','PW');
load('WorkingSignal.mat','signal_ifft','SIGNAL_ifft');
load('AnalogCircuitProcessing_1.mat','TVG_gain','Fixed_gain','FR','TT','FS');
%% Pulse compression frequency
FP = 2*BW;
TTT = TT(1:round(FS/FP):end);

%% Matching signal
SIGNAL_ifft = SIGNAL_ifft(1:round(FR/FS):end);
signal_ifft = signal_ifft(1:round(FR/FS):end);
Matching_signal = zeros(1,length(SIGNAL_ifft));
Matching_signal(1:length(signal_ifft)) = conj(fliplr(signal_ifft));
MatchingSignal_ABandPass = filter(AChebyshevTypeII_BandPass,Matching_signal);
TVG_Gain = zeros(1,length(MatchingSignal_ABandPass));
TVG_start = floor(0.5*length(TVG_gain));
TVG_Gain((1:length(signal_ifft))) = TVG_gain(TVG_start:TVG_start+length(signal_ifft)-1);
MatchingSignal_TVG = MatchingSignal_ABandPass.*TVG_Gain;
MatchingSignal_Fixed = MatchingSignal_TVG*Fixed_gain;
MatchingSignal_AD = MatchingSignal_Fixed;
MatchingSignal_DBandPass = filter(DChebyshevTypeII_BandPass,MatchingSignal_AD);
MatchingSignal_downB = MatchingSignal_DBandPass.*(exp(-1i*2*pi*FC*TT));
MatchingSignal_DLowPass = filter(DChebyshevTypeII_LowPass,MatchingSignal_downB);
MatchingSignal_downSB = MatchingSignal_DLowPass(1:round(FS/FP):end);
MatchingSignal = MatchingSignal_downSB/max(MatchingSignal_downSB);

% % % signal = real(MatchingSignal);
% % % Signal_fft = 2*fft(signal);Signal__FFT = abs(Signal_fft/length(signal));
% % % FD = Signal__FFT(1:length(signal)/2+1);FD(2:end-1) = FD(2:end-1);FD_F = FS*(0:(length(signal)/2))/length(signal);
% % % figure(1);scrsz = [20,40,1500,700];set(gcf,'Position',scrsz);
% % % subplot(2,1,1);plot(signal);set(gca,'FontSize',14);ylim([min(signal),max(signal)]);
% % % subplot(2,1,2);plot(FD_F,FD);set(gca,'FontSize',14);

%% Digital signal processing
for transmitter = TRANSMITTER
% transmitter = 222; 
    
    % 0) Digital signal
    filename = ['AnalogCircuitProcessing_',num2str(transmitter),'.mat'];
    load(filename,'Asignal','TT','FS');
%     DSignal = sum(Asignal);
    DSignal = Asignal(2,:);
    
    % 1) Bandpass filtering
    Dsignal_BandPassFilter = filter(DChebyshevTypeII_BandPass,DSignal); 
    % 2) Drop baseband        
    Signal_downB = Dsignal_BandPassFilter.*(exp(-1i*2*pi*FC*TT));
    % 3) Low pass filtering
    Signal_LowPass = filter(DChebyshevTypeII_LowPass,Signal_downB); 
    % 4) Bandpass down sampling
    Signal_downSB = Signal_LowPass(1:round(FS/FP):end);
    % 5) Matchedfilter
    Signal_matching = conv(Signal_downSB,MatchingSignal);
	Zero_Offset = floor(PW*FP)+1;  
    DSignal_Matching = Signal_matching(Zero_Offset:Zero_Offset+length(Signal_downSB)-1);  
	
    DATA = abs(DSignal_Matching);
    
    % Save data
    filename = ['SBSSS_',num2str(transmitter),'.mat'];
    save(filename,'MatchingSignal','DSignal','Dsignal_BandPassFilter',...
        'Signal_downB','Signal_LowPass','Signal_downSB',...
        'DATA','TTT','FP');
    disp(filename);
    sound(sin(2*pi*20*(1:4000)/100));
    
end

%% Save data
message = [num2str(transmitter),' groups of digital processing have been synthesized'];
disp(message);
sound(sin(2*pi*10*(1:4000)/100));

