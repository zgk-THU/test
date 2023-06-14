clear
close all
Nx=50;
Ny=50;
x=linspace(0,2*pi,Nx);
y=linspace(0,2*pi,Ny);
[Mx,My]=meshgrid(x,y);
Image_Plane_M=ones(2,2);
Image_Plane_Phase=2*pi*rand(2,2);
Image_Plane_Field=Image_Plane_M.*exp(1i*Image_Plane_Phase);

Fourier_Plane_Field=fft2(Image_Plane_Field);
Fourier_Plane_M=abs(Fourier_Plane_Field);
error=zeros(Nx,Ny);
for row=1:Ny
    for col=1:Nx
        error(row,col)=err(Mx(row,col),My(row,col),Image_Plane_Field,Fourier_Plane_M);
    end
end
% phase_show(Image_Plane_Field);
% imshow(Image_Plane_Field)
% title('phase angle (
% surf(Mx,My,error);
contourf(Mx,My,sqrt(error),40)
xlabel('x');
ylabel('y')
hold on
plot(Image_Plane_Phase(1,1),Image_Plane_Phase(1,2),'or');
colorbar
% colorbar('Ticks',[1e-2,pi/2,pi,3*pi/2,2*pi-1e-2],...
%          'TickLabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
% xlabel('x');ylabel('y');
% zlabel('error');
% Phase_guess=2*pi*rand(Nx,Ny);
% Image_guess=Image_Plane_M.*exp(1i*Phase_guess);
% img_it=Image_guess;
% for it=1:200
%     img_it=Proj_FM(img_it,Fourier_Plane_M);
%     img_it=Proj_IM(img_it,Image_Plane_M);
%     SE(it)=sum(sum((abs(fft2(img_it))-Fourier_Plane_M).^2));
% end
% plot(1:it_max,SE)
% % imshow(fftshift(abs(fft2(img_it))./max(max(abs(fft2(img_it))))));
function ER=err(x,y,IF,FM)
    XF=IF;
    XF(1,1)=exp(1i*x);
    XF(1,2)=exp(1i*y);
    FMM=abs(fft2(XF));
    ER=sum(sum(abs(FMM-FM).^2));
end
% % imshow(fftshift(M_show));
% image2=ifft2(ifftshift(image1));
% imshow(image1)
% for t=1:1000
%     imangle=angle(image2);
%     image=exp(1i*imangle);
%     image=fftshift(fft2(image));
%     imgabs=abs(image)./max(max(abs(image)));
%     corrcoef(x_true,imgabs);
%     if sim(1,2)>=0.9995
% 
%         break
%     else
%     imangle=angle(image2);
%     image2=exp(1i*imangle);
%     image3=fftshift(fft2(image2));
    % RMS_GS=zeros(500,1);
    % for k=1:512
    %     Phase_inc(:,k)=w(k);
    % end
    % AP_inc=Amp_inc.*exp(1i*Phase_inc);
    % AP_img=fftshift(fft2(AP_inc));
    % rec=ifft2(ifftshift(AP_img));
    % PH=mod(2*pi+angle(rec),2*pi);
    % imshow(PH./(2*pi))
    % Amplitude=abs(rec)./(max(max(abs(rec))));
    % imshow(abs(Amplitude)./(max(max(abs(Amplitude)))));
    % Amplitude=abs(AP_img)./(max(max(abs(AP_img))));
    % imshow(abs(Amplitude)./(max(max(abs(Amplitude)))));
    % % IFAR=4*A_1*(I_max/I_foot).^(2/5)
    % % % u_max=10^7*sqrt(0.7*A_1*alpha^0.6* ...
    % %     1^(4/15)*(I_max/I_foot)^0.4)
    % % E_L=4*pi*R_0^3*I_max/u_max
    % % rho_stag=16*rho_1*IFAR*(I_max/I_foot)^0.4
    % % amplication=P_stag/P_max
    % % rhoR_stag=rho_1*delta_1*IFAR^(2/3)*A_1^(2/3)*(I_max/I_foot)^(4/9)
    % % P_stag=P_max*IFAR^2 
    % % Ts = 1/50;
    % % t = 0:Ts:10-Ts;                     
    % % x = sin(2*pi*15*t) + sin(2*pi*20*t);
    % % 
    % % y=fft(x);
    % % fs=1/Ts;
    % % f = (0:length(y)-1)*fs/length(y);
    % % plot(f,abs(y))
    % % xlabel('Frequency (Hz)')
    % % ylabel('Magnitude')
    % % title('Magnitude')
    % clear
    % clc
    % t=linspace(-10,10,5000);
    % f=zeros(1,size(t,2));
    % omega0=10;
    % delta_omega=1;
    % N=10;
    % %phi=zeros(1,size(t,2));
    % % phi=2*pi*((1:N)./N).^2;
    % %phi=2*pi*((1:N)./N).^0.5;
    % phi=rand(1,size(t,2))*2*pi;
    % 
    % for k=1:N
    %     f=f+cos((omega0+(k-N/2)*delta_omega)*t+phi(k));
    % end
    % plot(t,f.^2)
    % 
    % 
