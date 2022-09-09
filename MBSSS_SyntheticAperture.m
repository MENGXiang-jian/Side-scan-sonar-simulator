%% Multi-subarray synthetic aperture focusing Synthetic Aperture
%% Document description

%% Initialization
clear;
close all;

c = 1500;
%% Load data
load('MSAF_ArrayManifold.mat','PROBE','Array_Position');
VerticalDistance = real(PROBE);
HeadingDistance = imag(PROBE);
load('ULADesigner.mat','ULA_SubarrayBeamAngle');
load('WorkingSignal.mat','DetectionDistance');

%% Synthetic aperture
SyntheticAperture = cell(1,length(ULA_SubarrayBeamAngle)-1); 
MAP = zeros(size(PROBE,1),size(PROBE,2));
SyntheticAperture_Complement = cell(1,length(ULA_SubarrayBeamAngle)-1); 
MAP_Complement = zeros(size(PROBE,1),size(PROBE,2)); 
for subarray = 1:length(ULA_SubarrayBeamAngle)-1
    
    map = zeros(size(PROBE,1),size(PROBE,2)); 
    % Geometry
    HeadingCore = Array_Position(1,subarray);
	Delt_Heading = abs(HeadingDistance-HeadingCore);
    ObliqueDistance = sqrt(Delt_Heading.^2+VerticalDistance.^2);
    HorizontalAngle = asind(Delt_Heading./ObliqueDistance);    
    % Logic
    map(ObliqueDistance <= DetectionDistance) = 1; 
    map(HorizontalAngle > (ULA_SubarrayBeamAngle(1,subarray)/2)) = 0;    
    
    SyntheticAperture{1,subarray} = map;
    MAP = MAP+map;    
    
    % Complement
    Complement = [135 102  68 35 1;...
                  167 134 101 67 34];
    map_Complement = map; 
%     map_Complement(Complement(1,subarray):Complement(2,subarray),1:size(PROBE,2)) = 1;    
    SyntheticAperture_Complement{1,subarray} = map_Complement;
    MAP_Complement = MAP_Complement+map_Complement;    

end
%% Save data
save('MSAF_SyntheticAperture.mat','SyntheticAperture','MAP','SyntheticAperture_Complement','MAP_Complement');
sound(sin(2*pi*10*(1:4000)/100));
%% Display MAP
figure(1)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

mesh(MAP)
view(2);
%% Display MAP_Complement
figure(2)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

mesh(MAP_Complement)
view(2);
%% Display Complement
figure(3)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

subarray = 3;
mesh(SyntheticAperture_Complement{1,subarray})
% mesh(Synthetic{1,subarray})
view(2);
%% Display MAP_Sea
Seabed_delt_xm = 0.01;    % 0.01
Seabed_delt_ym = 0.001;  % 0.0001
load('ULADesigner.mat','ULA_SubarrayLength');

[VerticalDistance,HeadingDistance] = meshgrid(0:Seabed_delt_xm:DetectionDistance,0:Seabed_delt_ym:sum(ULA_SubarrayLength(1:5)));
MAP_Sea = zeros(size(HeadingDistance,1),size(HeadingDistance,2));
for subarray = 1:length(ULA_SubarrayBeamAngle)-1
    
    map_Sea = zeros(size(HeadingDistance,1),size(HeadingDistance,2));
    
    % Geometry
    HeadingCore = Array_Position(1,subarray);    
    Delt_Heading = abs(HeadingDistance-HeadingCore);
    ObliqueDistance = sqrt(Delt_Heading.^2+VerticalDistance.^2);
    HorizontalAngle = asind(Delt_Heading./ObliqueDistance);      
    % Logic
    map_Sea(ObliqueDistance <= DetectionDistance) = 1; 
    map_Sea(HorizontalAngle > (ULA_SubarrayBeamAngle(1,subarray)/2)) = 0;
    
    MAP_Sea = MAP_Sea+map_Sea;
    
    disp(subarray);
end

save('MSAF_SyntheticAperture.mat','MAP_Sea','-append');
sound(sin(2*pi*10*(1:4000)/100));

figure(4)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

imagesc(MAP_Sea);
view(2);


