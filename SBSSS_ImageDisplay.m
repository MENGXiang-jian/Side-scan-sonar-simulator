%% SBSSS imaging
%% Tips


%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','TRANSMITTER','DR','Tracks');
load('ULADesigner.mat','ULA_SubarrayBeamAngle');

%% Imaging
Imaging_ratio = 0.5; % Single frame width ratio 0.65
ImageDelt_ym = 0.01;
Trajectory = round(Tracks(TRANSMITTER)/ImageDelt_ym)+1;
% Imaging_width = pi/180*DR*(ULA_SubarrayBeamAngle(1,1)/(length(ULA_SubarrayBeamAngle)-1));
Imaging_width = pi/180*DR*(ULA_SubarrayBeamAngle(1,2));
Imaging_width = round((Imaging_ratio*Imaging_width)/ImageDelt_ym);

for transmitter = TRANSMITTER

    Image_head = Trajectory(1,transmitter)-floor(Imaging_width/2);
    Image_tail = Image_head+Imaging_width-1;
    if Image_head <= 0
        Image_head = 1;
    end
    
    filename = ['SBSSS_',num2str(transmitter),'.mat'];
    load(filename,'DATA','TTT');
    Image = DATA;

    SingleFrameImaging = repmat(Image,(Image_tail-Image_head+1),1); 
    
    Image_imagesc(Image_head:Image_tail,:) = SingleFrameImaging;
    Image_pcolor(transmitter,:) = Image;    
    
    disp(transmitter);
end

%% Display
% (1) imagesc display
figure(1);
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);
X = [0 max(TTT)*c/2];
Y = [0 Tracks(max(TRANSMITTER))];
Image_imagesc_dB = 10*log10(Image_imagesc/max(max(Image_imagesc)));
% Image_imagesc_dB = Image_imagesc;
imagesc(X,Y,Image_imagesc_dB);
% shading interp;
colormap(hot);
caxis([-35 0]);

% (2) pcolor display
figure(2);
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);
X = [0 max(TTT)*c/2];
Y = [0 Tracks(max(TRANSMITTER))];
Image_pcolor_dB = 10*log10(Image_pcolor/max(max(Image_pcolor)));
% Image_pcolor_dB = Image_pcolor;
pcolor(TTT*c/2,Tracks(TRANSMITTER),Image_pcolor_dB);
xlabel('Transverse direction（m）','FontSize',24); 
ylabel('Navigation direction（m）','FontSize',24);
set(gca,'FontSize',24);
shading interp;
colormap(hot);
caxis([-28 0]);

%% Save data
save('SBSSS_Image.mat','Image_imagesc','Image_imagesc_dB','Image_pcolor','Image_imagesc_dB');
sound(sin(2*pi*10*(1:4000)/100));

