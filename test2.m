close all
d=1;
x=linspace(0,d,100);
qx=1;
y=sin(qx*pi/d*x);
plot(x,sin(1*pi/d*x).^2,x,sin(10*pi/d*x).^2)
