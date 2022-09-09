figure(1);

surf(Seabed_X(1,:),Seabed_Y(:,1),Seabed_Z);
shading interp;
xlabel('Transverse direction（m）','FontSize',16); 
ylabel('Navigation direction（m）','FontSize',16);
zlabel('Water depth（m）','FontSize',16);
set(gca,'ZDir','reverse','FontSize',16);
zlim([0,40]);  
% axis equal;


