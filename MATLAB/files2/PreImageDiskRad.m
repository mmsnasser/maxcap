function map = PreImageDiskRad (Lc,Lk,r,n,tol,Maxiter)
%
%
%
%%
t        = (0:2*pi/n:2*pi-2*pi/n).'; 
m        =  length(Lc)-1;
angk     =  angle(Lc);
cent     =  Lc;
radx     = (1-0.5*r).*Lk;
rady     =  r.*radx;
%%
thet(1:n,1)  =  pi/2;
for k=2:m+1
    thet(1+(k-1)*n:k*n,1) =  0;
end
%%

%%
alpha       =  0;
et(1:n,1)   =   exp(i.*t);
etp(1:n,1)  =   i.*exp(i.*t);
err = inf;
itr = 0;
while (err>tol)
    itr  =itr+1;  
    %
    for k=2:m+1
        et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*angk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
        etp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*angk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
    end
    %
    %
    A       =  exp(i.*(pi/2-thet)).*(et-alpha);
    for k=1:m+1
        gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*clog(et((k-1)*n+1:k*n,1)-alpha)); 
    end
    %%
    [mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
    h0            =  sum(h(1:n,1))/n;
    c             =  exp(-h0);
    R             =  h0.*sin(thet)-h;
    %
    fnet          = (gam+h+i.*mun)./A;
    zet           =  c.*(et-alpha).*exp((et-alpha).*fnet);
    rotwn         =  exp(-i.*R).*zet;
    for k=2:m+1
        Rk(k,1)     =  sum(R((k-1)*n+1:k*n,1))/n; 
        wnL         =  rotwn((k-1)*n+1:k*n,1);
        centk(k,1)  =  exp(i.*Rk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
        radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
    end
    cent   =  cent-1.0.*(centk-Lc);
    radx   =  radx-(1-0.5*r).*(radk-Lk) ;
    angk   =   angle(cent);
    rady   =  r.*radx;
    err    = (norm(centk-Lc,1)+norm(radk-Lk,1))/m;
    [itr err]
    error (itr,1) = err;
    itrk  (itr,1) = itr;
    %
figure(10);
clf; hold on; box on; axis equal
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-r','LineWidth',1.5);
    crv = zet(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
end 
drawnow
    %
    if itr>=Maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%%
for k=2:m+1
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*angk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*angk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
end
%
A       =  exp(i.*(pi/2-thet)).*(et-alpha);
for k=1:m+1
    gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*clog(et((k-1)*n+1:k*n,1)-alpha)); 
end
% 
[mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
h0            =  sum(h(1:n,1))/n;
c             =  exp(-h0);
fnet          = (gam+h+i.*mun)./A;
zet           =  c.*(et-alpha).*exp((et-alpha).*fnet);
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
map.cent  =  cent;
map.radx  =  radx;
map.c     =  c;
map.alpha =  alpha;
%
end