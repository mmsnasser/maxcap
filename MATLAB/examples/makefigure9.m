clear;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
global hr Rd cap et n ak
% 
%%
Rd    =  0.75;
%
rrv =  0.2:0.1:2;
%
hdisv = [];
akm   = [];
%% 
for jj=1:length(rrv)
hr  =  [rrv(jj); 0.2; 0.2; 0.2 ; 0.2 ; 0.2];
%
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
m    = length(hr);
x0c  = 0.75*exp([0:2/m:2-2/m]*pi*i);
x    = capmaxoptim(x0c);
% 
%
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
m    = length(hr);
z(1)  =  x(1);
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
hdisv = [hdisv hdis];
hdisv
akm   = [akm 2*pi*ak];
% 
%
[zet,zetp] = hyppolygonp(z,2^8);
%%
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
    else 
        plot(real(crv),imag(crv),'r-','LineWidth',1.5);
    end
end
%
plot(real(zet),imag(zet),':m','LineWidth',1.25);
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.25)
% plot(real(alpha),imag(alpha),'pr','LineWidth',1.5)
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
ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
%%
end
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
colr = ['r','b','k'];
for k=1:3
    plot(rrv,hdisv(k,:),colr(k),'LineWidth',1.5);
end
%
set(gca,'FontSize',14)
% axis equal
axis square
% axis([0.2 2 2.2 3])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
Leg=legend({'$\rho(z_1,z_2)$','$\rho(z_2,z_3)$','$\rho(z_3,z_4)$'},...
        'Interpreter','LaTeX','location','southwest');
xlabel('$r_1$','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigD6disv
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
box on; hold on;
plot(rrv,akm(1,:),'-ok','LineWidth',1.5);
plot(rrv,akm(2,:),'-vr','LineWidth',1.5);
plot(rrv,akm(3,:),'-^b','LineWidth',1.5);
plot(rrv,akm(4,:),'-+m','LineWidth',1.5);
%
set(gca,'FontSize',14)
% axis equal
axis square
% axis([0.2 2 2.2 3])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
Leg=legend({'$b_1$','$b_2$','$b_3$','$b_4$'},...
        'Interpreter','LaTeX','location','northwest');
xlabel('$r_1$','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigD6akv
%%