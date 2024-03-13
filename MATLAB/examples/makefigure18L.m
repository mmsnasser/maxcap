clear all;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
%
n       =   2^12;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
% 
hr    =  [1  ; 1  ; 1  ; 1  ; 1 ]; 
delt  =  [1  ; 1  ; 1  ; 1  ; 1 ];
deltv =  [0  ; delt ];
m     =  length(hr);
%
Rd    =  0.98;
xr    =  tanh(hr(1)/4)+0.02;
hzv   =  [linspace(xr,Rd,21)].';
%%
LB    =  5*(2*pi)./mu(tanh(hr(1)/2).^5);
UB    =  5*(2*pi)./mu(tanh(hr(1)/2));
%%
% 
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
% 
for kk=1:length(hzv)
    %
    for jj=1:m
        hz(jj,1) = hzv(kk)*exp(2i*pi*jj/m);
    end
    %
%
for k=1:m
    [Lc(k,1),crd(k,1)]   = HypDisk(hz(k),hr(k)/2);
end
Lk = 2*crd;
%
map    =  PreImageDiskRad ([0;Lc],[0;Lk],1,n,1e-12,200);
%% 
zet     =   map.zet;
zetp    =   map.zetp;
% 
et      =   map.et;
etp     =   map.etp;
% 
%%     
figure(11);
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on; axis equal
for k=1:m+1
    crv = zet(1+(k-1)*n:k*n);crv(n)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
end
%
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n)=crv(1);
    plot(real(crv),imag(crv),':r','LineWidth',1.5);
end
drawnow
%%
cent    =  map.cent; cent(1)=[];
alphav  =  cent;
alpha   =  0;
[cap,~] =  capm(et,etp,alphav,deltv,m,alpha);
cap
capv(kk,1) = cap;
%%
end
%% 
[hzv capv]
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(hzv,capv,'-b','LineWidth',1.5);
plot([0 1],[UB UB],':k','LineWidth',2);
plot([0 1],[LB LB],'--k','LineWidth',2);
set(gca,'FontSize',14)
axis square
axis([0.2 1 5.0 16])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
xlabel('$x$','Interpreter','latex')
ylabel('Capacity','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigS5v1
%%
