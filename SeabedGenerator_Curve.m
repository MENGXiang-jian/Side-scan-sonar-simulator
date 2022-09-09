%% Seabed model
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameter
load('START.mat','NavigationSpeed_mpers','Voyage');
%% Build model
Seabed_Xm = 260;
Seabed_Ym = Voyage;
Seabed_Zm = 30;

Seabed_delt_xm = 0.1;
Seabed_delt_ym = 0.01; 

Seabed_x = round(Seabed_Xm/Seabed_delt_xm)+1;
Seabed_y = round(Seabed_Ym/Seabed_delt_ym)+1;

Seabed_Z = Seabed_Zm*ones(Seabed_y,Seabed_x);
[Seabed_X,Seabed_Y] = meshgrid(0:Seabed_delt_xm:Seabed_Xm,0:Seabed_delt_ym:Seabed_Ym);

%% Seabed target
CURVE = real(30*exp(1j*2*pi*0.1*Seabed_Y(:,1)));
curve_X = [35 80 125 170 215];
for yy = 1:size(Seabed_Y,1)    
    for curve = 1:length(curve_X)
        xx = round(curve_X(curve)/Seabed_delt_xm)+1+round(CURVE(yy)/Seabed_delt_xm);
        zz = Seabed_Zm-2;
        Seabed_Z(yy,xx) = zz;  
    end
end
% % % % Random point
% % % num = 3000;aa = numel(Seabed_Z);
% % % Seabed_Z(datasample((find(Seabed_Z == Seabed_Zm))',num)) = Seabed_Zm-rand(1,num)*2;

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








