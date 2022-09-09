%% Multi-subarray synthetic aperture focusing array manifold
%% Document description

%% Initialization
clear;
close all;

c = 1500;
%% Single frame imaging
% Imaging size
load('WorkingSignal.mat','DetectionDistance');
Image_side = DetectionDistance;
Image_width = 1.67;
% Pixel
load('Seabed.mat','Seabed_delt_xm','Seabed_delt_ym');
ImageDelt_xm = 0.1;   Seabed_delt_xm;
ImageDelt_ym = 0.01;  Seabed_delt_ym;
Image_x = 0:ImageDelt_xm:Image_side;  Probe_x = 0.5*ImageDelt_xm:ImageDelt_xm:Image_side;
Image_y = 0:ImageDelt_ym:Image_width; Probe_y = 0.5*ImageDelt_ym:ImageDelt_ym:Image_width;

% IMAGE
[Image_X,Image_Y] = meshgrid(Image_x,Image_y);
IMAGE = Image_X+1j*Image_Y;
% PROBE
[Probe_X,Probe_Y] = meshgrid(Probe_x,Probe_y);
PROBE = Probe_X+1j*Probe_Y;
%% Array manifold matrix
load('ULADesigner.mat','ULA_SubarrayPosition');
Array_Position = ULA_SubarrayPosition(1:5)+0.5*Image_width;
load('START.mat','POSITION','SUBARRAY');

TOA = cell(size(PROBE,1),size(PROBE,2));
for position = POSITION 
%     position = 10;
    
    % Universal        
    for probe = 1:numel(PROBE)
        PROBE_X = real(PROBE(probe));
        PROBE_Y = imag(PROBE(probe));
        toa = zeros(1,length(ULA_SubarrayPosition)-1);
        % Driving vector        
        for subarray = SUBARRAY
            toa(1,subarray) = ((sqrt((PROBE_Y-Array_Position(1,2))^2+PROBE_X^2))/c)+(sqrt((PROBE_Y-Array_Position(1,subarray))^2+PROBE_X^2))/c;                  
        end
        TOA{probe} = toa;
    end
    
    % Replace
    TOA_replace = TOA;
    filename = ['NavigationScanning_',num2str(position),'.mat'];
    load(filename,'Photo');
    photo = Photo(Photo ~= 0);
    for pp = 1:length(photo)              
        photo_y = floor(real(photo(pp))/ImageDelt_ym)+1;
        photo_x = floor(imag(photo(pp))/ImageDelt_xm)+1;       
        
        toa_replace = zeros(1,length(ULA_SubarrayPosition)-1);
        % Driving vector        
        for subarray = SUBARRAY
            toa_replace(1,subarray) = ((sqrt((imag(photo(pp))-Array_Position(1,2))^2+real(photo(pp))^2))/c)+(sqrt((imag(photo(pp))-Array_Position(1,subarray))^2+real(photo(pp))^2))/c;                  
        end
        TOA_replace{photo_y,photo_x} = toa_replace; % Take only one group
        
        message = ['位置',num2str(position),'完成替换处理'];
        disp(message);
    end    
    
    % Save data       
    filename = ['MSAF_ArrayManifold_',num2str(position),'.mat'];
    save(filename,'TOA','TOA_replace');
    sound(sin(2*pi*20*(1:4000)/100));
    disp(position);
end

%% Save data
save('MSAF_ArrayManifold.mat','ImageDelt_xm','ImageDelt_ym','Probe_x','Probe_y','Image_x','Image_y','IMAGE','PROBE','Array_Position');
sound(sin(2*pi*10*(1:4000)/100));
%% Display


