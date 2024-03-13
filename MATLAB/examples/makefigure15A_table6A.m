clear all;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
global hr Rd cap et n fig1 videow
% 
%%
Rd    =  0.75;
% 
hr  =  [0.4; 0.4; 0.4; 0.2 ; 0.2];
%
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
%%
% 
% 
% videow = VideoWriter('maxD4','MPEG-4'); % Video ========
% videow.Quality = 99; % Video ========
% videow.FrameRate = 5; % Video ========
% open(videow); % Video ========
fig1=figure(1);
%
%
%
m    = length(hr);
x0c  = [-0.7 -0.3  0.2  0.5  0.7];
x = capmaxoptimre(x0c);
% 
%%
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
m    =  length(hr);
z    =  x;
% 
[as,bs] = sort(z);
zv = z(bs); 
hrv = hr(bs); 
%
hdis = [];  cnthd = []; 
for k=1:m-1
    hdis   = [hdis  ; rho(zv(k),zv(k+1))];
    cnthd  = [cnthd ; rho(zv(k),zv(k+1))-hrv(k)-hrv(k+1)];
end
%
format short g
z.'
[hdis cnthd]
% 
%        1.2558      0.45576
%        1.2719      0.47188
%       0.86372      0.26372
%       0.50046      0.10046
%%
%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5);
for k=2:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    if hr(k-1)==0.2
        plot(real(crv),imag(crv),'b-','LineWidth',1.5);
    elseif hr(k-1)==0.4
        plot(real(crv),imag(crv),'r-','LineWidth',1.5);
    end
end
%
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.0)
% plot(real(alpha),imag(alpha),'pk','LineWidth',1.5)
str = sprintf('cap$=$ %1.6f ', -cap);
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
print -depsc FigReD6n31
%%
