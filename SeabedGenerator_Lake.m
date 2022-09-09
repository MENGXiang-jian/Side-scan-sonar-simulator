%% Seabed model
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameter
load('START.mat','NavigationSpeed_mpers','Voyage');
%% Build model
Seabed_Xm = 100;
Seabed_Ym = Voyage;
Seabed_Zm = 30;

Seabed_delt_xm = 0.1;
Seabed_delt_ym = 0.01; 

Seabed_x = round(Seabed_Xm/Seabed_delt_xm)+1;
Seabed_y = round(Seabed_Ym/Seabed_delt_ym)+1;

Seabed_Z = Seabed_Zm*ones(Seabed_y,Seabed_x);
[Seabed_X,Seabed_Y] = meshgrid(0:Seabed_delt_xm:Seabed_Xm,0:Seabed_delt_ym:Seabed_Ym);

%% Lake
SC_Ym = [5 10 15];
SC_Xm = [20];

for yy = 1:length(SC_Ym)
    for xx = 1:length(SC_Xm)        
        % target position
        SC_ym = SC_Ym(yy);
        SC_xm = SC_Xm(xx);
        % target size
        Y_SC = 4;
        X_SC = 2;
        SC_Zm = Seabed_Zm-10;        
        % Place target
        SC_xs = round((SC_xm-Y_SC)/Seabed_delt_xm)+1;
        SC_xe = round((SC_xm+Y_SC)/Seabed_delt_xm)+1; 
        SC_ys = round((SC_ym-0.5*X_SC)/Seabed_delt_ym)+1; 
        SC_ye = round((SC_ym+0.5*X_SC)/Seabed_delt_ym)+1;   
        
        Seabed_Z(SC_ys:SC_ye,SC_xs:SC_xe) = SC_Zm;     
    end
end


% Line
Line_X = 50:10:90;
for line_x = 1:length(Line_X)
    Seabed_Z(:,Line_X(line_x)/Seabed_delt_xm+1:Line_X(line_x)/Seabed_delt_xm+1+10) = Seabed_Zm-1;
end
Line_Y = 2:2:18;
for line_y = 1:length(Line_Y)
    Seabed_Z(Line_Y(line_y)/Seabed_delt_ym+1:Line_Y(line_y)/Seabed_delt_ym+1+10,50/Seabed_delt_xm+1:end) = Seabed_Zm-1;
end

%% Save data
save('Seabed.mat','Seabed_Zm','Seabed_X','Seabed_Y','Seabed_Z','Seabed_delt_xm','Seabed_delt_ym');
sound(sin(2*pi*10*(1:4000)/100));

%% Display
figure(1);
% scrsz = [20,40,1500,700];
% set(gcf,'Position',scrsz);

surf(Seabed_X(1,:),Seabed_Y(:,1),Seabed_Z);
shading interp;
xlabel('Transverse direction（m）','FontSize',15); 
ylabel('Navigation direction（m）','FontSize',15);
zlabel('Water depth（m）','FontSize',15);
set(gca,'ZDir','reverse','FontSize',15);
set(gca,'XTick',(0:50:255));
set(gca,'YTick',(0:4:20));
zlim([15,30]);  
% axis equal;








