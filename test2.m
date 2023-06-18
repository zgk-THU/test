close all
clc
d=1;
x=linspace(0,d,100);
y=linspace(0,d,100);
dx=x(2)-x(1);
dy=y(2)-y(1);
[X_mesh,Y_mesh]=meshgrid(x,y);
energy=@(U) sum(sum(U.^2))*dx*dy;
mode=@(qx,qy) (2/d)*sin(qx*pi/d*X_mesh).*...
    sin(qy*pi/d*Y_mesh);
% energy(mode(2,2)+mode(1,1))
Uxy=2*mode(1,1)+1*mode(2,2)+mode(3,2);
figure
set(gcf,'Position',[1000,100,500,500])
t=tiledlayout(1,1);
nexttile
t.Padding='compact';
imagesc(Uxy)
colormap twilight
% colormap copper(128)
colorbar
