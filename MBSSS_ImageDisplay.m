%% MBSSS ImageDisplay
%% Document description

%% Initialization
clear;
close all;

c = 1500;
%% Load data
load('MSAF_ArrayManifold.mat','ImageDelt_xm','ImageDelt_ym');

%% Single frame imaging
load('MSAF_ArrayManifold.mat','Image_x','Image_y','IMAGE','PROBE','Array_Position');

%% Navigation information (No correction)
load('NavigationTrajectory.mat','Trajectory');
Trajectory_pixel = round(Trajectory/ImageDelt_ym)+1;

%% MBSSS imaging
load('START.mat','POSITION');
for position = POSITION
%     position = 10;

    filename = ['MBSSS_Beamforming_',num2str(position),'.mat'];
    load(filename,'DATA');
    
    Image_head = Trajectory_pixel(1,position)-floor(size(DATA,1)/2);
    Image_tail = Image_head+size(DATA,1)-1;
    if Image_head <= 0
        Image_head = 1;
    end 
    
    % Signal frame imaging
    SingleFrameImaging = DATA(end-(Image_tail-Image_head):end,:);
    
    % Mosaic imaging
    Image(Image_head:Image_tail,:) = SingleFrameImaging;    
    
    disp(position);
end
Image = Image.^2;
%% Display
% (1) data display
figure(1)
scrsz = [20,40,600,340];
set(gcf,'Position',scrsz);
X = [0 250];
Y = [0 20];
Image_imagesc_dB = 10*log10(Image/max(max(Image)));
imagesc(X,Y,Image_imagesc_dB);
xlabel('Transverse distance（m）','FontSize',10); 
ylabel('Navigation distance（m）','FontSize',10);
set(gca,'YDir','normal'); 
set(gca,'FontSize',10);
shading interp;
colormap(hot);
caxis([-38 0]); 

% (2) pcolor display
figure(2)
scrsz = [20,40,1500,600];
set(gcf,'Position',scrsz);
Image_imagesc_dB = 10*log10(Image/max(max(Image)));
pcolor(Image_imagesc_dB);
xlabel('Transverse distance（m）','FontSize',15); 
ylabel('Navigation distance（m）','FontSize',15);
set(gca,'FontSize',24);
shading interp;
colormap(hot);
caxis([-38 0]); 

%% Save data
save('MBSSS_Image.mat','Image','Image_imagesc_dB');
sound(sin(2*pi*10*(1:4000)/100));

%% Display


