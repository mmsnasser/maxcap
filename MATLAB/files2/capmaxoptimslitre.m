function x = capmaxoptimslitre(x0c)
% 
%
%%
% global hr Rd %cap et zet n fig1 
% 
%%
% 
%
%
%%
% 
% 
%
%
%
fun     = @capfunmaxslitre;
nonlcon = @mindistance;
%
x0   = x0c;
%
% 
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = -1+zeros(size(x0));%-inf;%[];
ub =  1+zeros(size(x0));% inf;%[];
options = optimoptions('fmincon','Display','iter','FunctionTolerance',1e-7,'Algorithm','interior-point','MaxIterations',200,'MaxFunctionEvaluations',400);
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
%
% 
%%
end


%%
function [c,ceq] = mindistance(x)
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
dis  =  0.05;
global hr Rd 
%
m   = length(hr);
z   = x;
ind = 0;
for k=1:m
    for j=1:k-1
        ind  = ind+1;
        c(ind) = hr(k)/2+hr(j)/2+dis-rho(z(k),z(j));
    end
end
for k=1:m
    ind = ind+1;
    c(ind) = abs(z(k))-Rd;
end
ceq = [];
end