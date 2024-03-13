clear all;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
global hr Rd cap et zet n  
% 
%%
Rd    =  0.75;
hr    = [0.4; 0.4; 0.4; 0.4; 0.4]; 
%
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
% 
%
%%
figure(1);
m    = length(hr);
x0c  = 0.65*linspace(-1,1,m);
x = capmaxoptimslitre(x0c);
%%
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
z    =  x;
%% 
[as,bs] = sort(z);
zv = z(bs); 
hrv = hr(bs); 
%
hdis = []; 
for k=1:m-1
    hdis   = [hdis  ; rho(zv(k),zv(k+1))];
end
%
'---------------------------'
sprintf('%1.6f\n', cap)
sprintf('%1.6f\n', z)
sprintf('%1.6f\n', hdis)
% ans =
%     '6.701089
%      '
% ans =
%     '-0.750000
%      -0.468629
%      0.000000
%      0.468629
%      0.750000
%      '
% ans =
%     '0.929285
%      1.016625
%      1.016624
%      0.929285
%      '
%      '
%%
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
        plot(real(crv),imag(crv),'b-','LineWidth',1.5);
    else 
        plot(real(crv),imag(crv),'r-','LineWidth',1.5);
    end
end
%
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
print -depsc FigSlit5e
%%