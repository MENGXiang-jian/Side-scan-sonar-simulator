%% Multi subarray synthetic aperture dynamic focusing beamforming
%% Document description

%% Initialization
clear;
close all;

c = 1500;
%% Load data
load('START.mat','F_Interpolated','F_Original','F_Sampling','POSITION');
load('WorkingSignal.mat','FC','BW','PW','DetectionDistance','DetectionPeriod','Tsignal');
load('ULADesigner.mat','ULA_SubarrayLength');
load('MSAF_SyntheticAperture.mat','SyntheticAperture','SyntheticAperture_Complement');
load('MBSSS_AnalogCircuitProcessing.mat','TVG_gain','Fixed_gain','ArtificialGain');

%% Matching signal
    Tsignal_Original = Tsignal(1:F_Interpolated/F_Original:end);
    TSignal_conj = conj(fliplr(Tsignal_Original));
    % Bandpass
    TSignal_BandPassFilter = filter(ABandPass_ChebyII,TSignal_conj);
    % TVG
    TVG_start = 100000;
    TSignal_TVG  = TSignal_BandPassFilter.*TVG_gain(TVG_start:TVG_start+length(TSignal_BandPassFilter)-1);        
    % Fixed gain
    TSignal_Fixed = TSignal_TVG*Fixed_gain/ArtificialGain;
    % AD
    TSignal_AD = TSignal_Fixed(1:F_Original/F_Sampling:end);
    % Bandpass filter
    Tsignal_BandPassFilter = filter(DBandPass_ChebyII,TSignal_AD); 
    % Drop baseband
    TSignal_downB = Tsignal_BandPassFilter.*(exp(-1i*2*pi*FC*(0:1/(F_Sampling*FC):(length(Tsignal_BandPassFilter)-1)*(1/(F_Sampling*FC)))));    
    % Low pass filtering
    TSignal_LowPass = filter(DLowPass_ChebyII,TSignal_downB);     
    % Bandpass down sampling
    TSignal_downSB = TSignal_LowPass(1:((F_Sampling*FC)/(BW*4)):end);
    Matchingsignal = zeros(1,2*length(TSignal_downSB));
    Matchingsignal(1,1:length(TSignal_downSB)) = TSignal_downSB;
    
%% Beamforming
tic
for position = POSITION
%     position = 10; 
    % Signal data
    filename = ['MBSSS_AnalogCircuitProcessing_',num2str(position),'.mat'];
    load(filename,'Dsignal','T_Dsignal');
    SIGNAL = Dsignal;
    % ArrayManifold data
    filename = ['MSAF_ArrayManifold_',num2str(position),'.mat'];
    load(filename,'TOA','TOA_replace');
    Toa = TOA;
    
    DATA = zeros(size(Toa,1),size(Toa,2));
    for probe = 1:numel(Toa)
        message = [num2str(position),':',num2str(probe),'/',num2str(numel(TOA)),' is probed'];
        disp(message);
        
% 1) Signal fragment        
        toa = Toa{probe};
        Toa_delta = toa-toa(1,3);
        % head
        toa_head = round(toa*F_Sampling*FC)+1;
%         toa_head = toa_head-round((0.5*PW)*F_Sampling*FC);  
        toa_head(toa_head <= 0) = 1;
        % tail
        toa_tail = (round(toa*F_Sampling*FC)+1)+round(PW*F_Sampling*FC); 
        toa_tail(toa_tail >= (round((2*DetectionDistance/c)*F_Sampling*FC)+1)) = round((2*DetectionDistance/c)*F_Sampling*FC)+1;
        % Intercepted fragment
        Toa_head = min(toa_head);    
        Toa_tail = max(toa_tail);        
        % Effective fragment    
        Toa_head_bias = round(toa(1,3)*F_Sampling*FC)+1-Toa_head+1;
        Toa_tail_bias = Toa_head_bias+PW*F_Sampling*FC;       
        bias = Toa_head_bias:Toa_tail_bias;        
% 2) Frequency domain beamforming FFT
        BFsignal_sum = zeros(1,Toa_tail-Toa_head+1);    
        aperture = zeros(1,length(ULA_SubarrayLength)-1);
        for subarray = 1:length(ULA_SubarrayLength)-1
            
            Synthetic_Aperture = SyntheticAperture_Complement{1,subarray}; % _Complement
            aperture(1,subarray) = Synthetic_Aperture(probe);
            if (Synthetic_Aperture(probe) == 1) && (Toa_head < Toa_tail)
                % (1) Signal fragment 
                signal = SIGNAL(subarray,Toa_head:Toa_tail);
%                 signal = a;
%                 plot(signal);
                % (2) FFT processing 
                N = size(signal,2); % Even number recommended
                FFTRsignal = 2*fft(signal,N);
                fftsignal = abs(FFTRsignal(1,1:floor(N/2)+1));
                f = F_Sampling*FC*(0:floor(N/2))/N;                
%                 plot(f,fftsignal);
                % (3) FFT point
                Narrowband = 0.01*FC;                                    % Minimum standards
                Narrowbandpoint = floor(Narrowband/(F_Sampling*FC)*N); % Span points
                NarrowbandPoint = 1;                                     % Real span points
                FFT_FL = 0*10^3; % Low frequency 390*10^3    (min = 0*10^3)
                FFT_FH = max(f); % High frequency 410*10^3 (max = max(f))
                FFT_FLpoint = round(FFT_FL/(F_Sampling*FC)*N)+1; % Low frequency point
                FFT_FHpoint = round(FFT_FH/(F_Sampling*FC)*N)+1; % High frequency point
                % Extraction frequency and frequency point
                FFTPoint = FFT_FLpoint:NarrowbandPoint:FFT_FHpoint;
                FFTfre = (FFTPoint-1)*(F_Sampling*FC)/N; 
                % (4) Phase
                toa_delta = Toa_delta(1,subarray);
                Signal = zeros(1,length(FFTfre));
                for fre = 1:length(FFTfre)            
                    phase = exp(1j*2*pi*toa_delta*FFTfre(fre)); 
                    Signal(1,fre) = FFTRsignal(1,FFTPoint(fre))*phase;
                    a = 1;
                end               
                BFsignal = zeros(1,N);
                BFsignal(1,FFTPoint) = Signal;          
                % (5) Sum
                BFsignal_sum = BFsignal_sum+BFsignal;  
            end
        end
% 3) IFFT        
        BFsignal= real(ifft(BFsignal_sum)); 
% 4) Effective fragment
        ShareEqually = sum(aperture);
        if ShareEqually == 0
            ShareEqually = 1;
        end       
        BFSignal = zeros(1,Toa_tail_bias);
        BFSignal(1,1:length(BFsignal)) = BFsignal/ShareEqually;       
        BFSIGNAL = BFSignal(1,Toa_head_bias:Toa_tail_bias);         
% 5) Bandpass filter
        Dsignal_BandPassFilter = filter(DBandPass_ChebyII,BFSIGNAL);
% 6) Drop baseband
        Signal_downB = Dsignal_BandPassFilter.*(exp(-1i*2*pi*FC*(0:1/(F_Sampling*FC):(length(Dsignal_BandPassFilter)-1)*(1/(F_Sampling*FC)))));
% 7) Low pass filtering
        Signal_LowPass = filter(DLowPass_ChebyII,Signal_downB);        
% 8) Bandpass down sampling
        Signal_downSB = Signal_LowPass(1:((F_Sampling*FC)/(BW*4)):end);          
% 9) Matchedfilter    
        MatchingSignal = Matchingsignal(1,1:length(Signal_downSB));
        Signal_matching = conv(Signal_downSB,MatchingSignal);
        Zero_Offset = floor(PW*BW*4)+1;
        DSignal_Matching = Signal_matching(Zero_Offset:Zero_Offset+length(Signal_downSB)-1);
% 10) Statistics
        data_fragment = round(5*0.5*(1/BW)*(4*BW));
        if data_fragment > length(DSignal_Matching)-1
            data_fragment = length(DSignal_Matching)-1;
        end
        DATA_signal = DSignal_Matching(1,1:data_fragment+1);
        DATA(probe) = sum(abs(DATA_signal)); 
        
%         test = 1;
    end
    % Save data
    filename = ['MBSSS_Beamforming_',num2str(position),'.mat'];
    save(filename,'DATA');
    sound(sin(2*pi*20*(1:4000)/100));
    disp(filename);
    toc
end
toc   
%% Save data

sound(sin(2*pi*10*(1:4000)/100));
%% Display


