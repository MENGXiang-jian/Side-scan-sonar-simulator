%% Waveform display
%% Initialization
% clear;
% close all;

c = 1500;
%% Load data
load('START.mat','FC','FS');
load('ArraySignal_7.mat','Signal');
load('AnalogCircuitProcessing_111.mat','Asignal_Original','Asignal_Transducer','Asignal_BandPassFilter',...
    'Asignal_TVG','Asignal','T','FR','TT','Noise','TVG_gain','Fixed_gain');
load('SBSSS_111.mat','MatchingSignal','DSignal','Dsignal_BandPassFilter','Signal_downB','Signal_LowPass',...
    'Signal_downSB','DATA','TTT','FP');
signal1 = Dsignal_BandPassFilter; % Dsignal(3,:)

load('AnalogCircuitProcessing_222.mat','Asignal_Original','Asignal_Transducer','Asignal_BandPassFilter',...
    'Asignal_TVG','Asignal','T','FR','TT','Noise','TVG_gain','Fixed_gain');
load('SBSSS_222.mat','MatchingSignal','DSignal','Dsignal_BandPassFilter','Signal_downB','Signal_LowPass',...
    'Signal_downSB','DATA','TTT','FP');
signal2 = Dsignal_BandPassFilter;

%% Display
% (1) Time domain analysis
figure(1);
plot(T,signal1);

hold on;

plot(T,signal2);

xlim([0,max(T)]); 
legend('Noise','Echo signal');
xlabel('Detection distance(m)','FontSize',15); 
ylabel('Amplitude(V)','FontSize',15); 
set(gca,'FontSize',15);

% (2) Frequency domain analysis
figure(2)
N = length(signal1);
FFTRsignal = 2*fft(signal1,N);
fftsignal = abs(FFTRsignal(1,1:floor(N/2)+1));
f = FS*(0:floor(N/2))/N; 
plot(f,fftsignal); 

hold on;

N = length(signal2);
FFTRsignal = 2*fft(signal2,N);
fftsignal = abs(FFTRsignal(1,1:floor(N/2)+1));
f = FS*(0:floor(N/2))/N; 
plot(f,fftsignal); 

legend('Noise','Echo signal');
xlabel('Frequence(Hz)','FontSize',15); 
ylabel('Amplitude(V)','FontSize',15); 
set(gca,'FontSize',15);
xlim([3*10^5,5*10^5]);

%% Save data

sound(sin(2*pi*10*(1:4000)/100));
%% Display


