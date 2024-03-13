clear all;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
global hr Rd cap et zet n fig1 
% 
%%
Rd    =  0.75;
hr    =  [1.6; 0.4; 0.4; 0.4; 0.4; 0.4]; 
%
n       =   2^12;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
% 
%
%%
fig1=figure(1);
m    = length(hr);
rndv = 0.75+0.25*rand(1,6);
x0c  = 0.6*rndv.*exp([0:2/m:2-2/m]*pi*i);
x = capmaxoptimrad(x0c);
%%
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
z(1) =  x(1);
for k=2:m
    z(k)  =  x(2*k-2)+i*x(2*k-1);
end
% 
thet = Arg(z,0);
[as,bs] = sort(thet);
zv = z(bs); zv(end+1)=zv(1);
hrv = hr(bs); hrv(end+1)=hrv(1);
%
hdis = [];   
for k=1:m
    hdis   = [hdis  ; rho(zv(k),zv(k+1))];
end
% 
format short g
hdis
%        2.8236
%        2.5089
%        2.4931
%        2.4931
%        2.5089
%        2.8236
%%
[xt,xtp] = hyppolygonp(z,2^8);
%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = zet((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5);
for k=2:m+1
    crv = zet((k-1)*n+1:k*n); crv(n)=crv(1);
    if hr(k-1)==0.4
        plot(real(crv),imag(crv),'b-','LineWidth',2);
    else 
        plot(real(crv),imag(crv),'r-','LineWidth',2);
    end
end
%
plot(real(xt),imag(xt),':m','LineWidth',1.0);
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.0)
% plot(real(alpha),imag(alpha),'pk','LineWidth',1.5)
str = sprintf('cap$=$ %1.6f ', cap);
title(str)
set(gca,'FontSize',14)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigRad6n1L
%%