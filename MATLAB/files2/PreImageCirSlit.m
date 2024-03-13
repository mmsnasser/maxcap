function map = PreImageCirSlit(rad,alp,bet,r,n,tol,Maxiter)
%
%
%
%%
t        = (0:2*pi/n:2*pi-2*pi/n).'; 
m        =  length(rad);
%%
rad    =  rad(:); alp = alp(:); bet = bet(:);
alpp   =  rad.*exp(i*alp);
betp   =  rad.*exp(i*bet);
seccp  = (alpp+betp)/2;
secL   =  abs(alpp-betp);
%
for k=1:m
    if bet(k)<alp(k)
        error();
    end
end
%
dskc   =  seccp;
dskr   =  abs(alpp-betp)/4;
%
%
%%
for k=1:m
    while abs(dskc(k))+dskr(k)>0.95
        dskr(k) = 0.9*dskr(k);
        if dskr(k)<1e-6
            error('===r is too small')
        end
    end
end    
%%
et(1:n,1)   =    exp(i*t);
etp(1:n,1)  =  i*exp(i*t);
%
err = inf;
itr = 0;
while (err>tol)
    itr  =itr+1;  
    %
    for k=1:m
        Jk = 1+k*n:(k+1)*n;
        et (Jk,1)   =  dskc(k)+dskr(k)*exp(-i*t);
        etp(Jk,1)   =       -i*dskr(k)*exp(-i*t);
    end
    %
    %
    A       =   et;
    gam     =  -log(abs(et));
    %
    [mun,h] = fbie(et,etp,A,gam,n,5,[],1e-14,100);
    %
    fnet    = (gam+h+i.*mun)./A;
    c       =  exp(-mean(h(1:n)));
    zet     =  c*et.*exp(et.*fnet);
    %
figure(1);
clf
hold on;box on; axis equal
crv = et(1:n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5);
plot(0,0,'pm','LineWidth',1.5);
%
plot(real(dskc),imag(dskc),'om','LineWidth',2);
plot(real(seccp),imag(seccp),'oc','LineWidth',2);
for k=1:m
    tt=linspace(alp(k),bet(k),50);
    crv = rad(k)*exp(i*tt);
    plot(real(crv),imag(crv),'k','LineWidth',2);
end
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    plot(real(zet(Jk)),imag(zet(Jk)),'r','LineWidth',1.0);
end
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    plot(real(et(Jk)),imag(et(Jk)),'b','LineWidth',1.0);
end    
drawnow 
%
    %    
    for k=1:m
        Jk         =  1+k*n:(k+1)*n;
        angk       =  carg(zet(Jk));
        alpk(k,1)  =  min(angk);
        betk(k,1)  =  max(angk);
        radk(k,1)  =  exp(mean(h(Jk))-mean(h(1:n)));
    end
    alpkp   =  radk.*exp(i*alpk);
    betkp   =  radk.*exp(i*betk);
    dirk    =  exp(i*(alpk+betk)/2);
    secLk   =  abs(alpkp-betkp);
    %
    %
    dskc  =  dskc + r.*(rad-radk).*dirk;   
    dskr  =  dskr + r.*((secL-secLk)/4);
    %
    for k=1:m
        while abs(dskc(k))+dskr(k)>0.999
        dskr(k) = 0.99*dskr(k);
        dskc(k) = 0.99*dskc(k);
        if dskr(k)<1e-6
            error('===r is too small')
        end
        end
    end
    %    
    err   =  (norm(alpp-alpkp,1)+norm(betp-betkp,1)+norm(rad-radk,1))/m;
    %
    %
    [itr err]
    errv (itr,1) = err;
    itrk  (itr,1) = itr;
    %
    if itr>=Maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%%
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    et (Jk,1)   =  dskc(k)+dskr(k)*exp(-i*t);
    etp(Jk,1)   =       -i*dskr(k)*exp(-i*t);
end
%
A       =   et;
gam     =  -log(abs(et));
%
[mun,h] = fbie(et,etp,A,gam,n,5,[],1e-14,100);
%
fnet    = (gam+h+i.*mun)./A;
c       =  exp(-mean(h(1:n)));
zet     =  c*et.*exp(et.*fnet);
%
for k=1:m+1
    Jk=(k-1)*n+1:k*n;
    zetp(Jk,1) = derfft(real(zet(Jk)))+i*derfft(imag(zet(Jk)));
end
%%
map.fnet  =  fnet;
map.zet   =  zet;
map.zetp  =  zetp;
map.et    =  et;
map.etp   =  etp;
map.cent  =  dskc;
map.rad   =  dskr;
map.c     =  c;
%%
end