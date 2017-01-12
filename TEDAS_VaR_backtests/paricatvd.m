%{
Function to estimate the final parameters for ICATVD using the 
GARCH assumption about the independent factors and 
the Normal Inverse Gaussian distribution 
%}
function [FuncVal,M2,M3,M4,skkurt,EstPar,RhoMat,ZetaMat] = paricatvd(InitPar,SigMat,F,AA,nflag)

if nargin==6
    if ~ischar(nflag)
        error('The input "flag" must be the character string ''forecast''.')
    else
        isforec = strcmpi(nflag,'forecast');
        if ~isforec
            error('Unrecognized "flag" input.')
        end
    end
else
    isforec = false;
end

F      = F(2:end,:);
[T,k]  = size(F);

HtLast  = diag(SigMat(T,:));
HtForec = diag(SigMat(T+1,:));

Alpha    = zeros(1,T+1);
Beta     = zeros(1,T+1);
Delta    = zeros(1,T+1);
Rho      = zeros(1,T+1);
Zeta     = zeros(1,T+1);

EstPar   = zeros(k,8);
AlphaMat = zeros(T,k);
BetaMat  = zeros(T,k);
DeltaMat = zeros(T,k);
ZetaMat  = zeros(T,k);
RhoMat   = zeros(T,k);
FVal     = zeros(1,k);

parinit  = [0.1 0.1 -0.1 1 -2 0.5 -0.9 0.01];

options = optimoptions('fminunc','MaxFunEvals',100000, 'MaxIter', 100000);

parfor j = 1:k
    fvec = F(:,j);
    fvol = SigMat(:,j);
    [estpar,fval] = fminunc(@(param)secstepLL(fvec,param,InitPar(j,:),fvol),parinit,options);
    EstPar(j,:)   = estpar;
    FVal(j)       = -fval;
    
    [~,rho,zeta,alpha,beta,~,~] = secstepLL(fvec,estpar,InitPar,fvol);
    
    Alpha = alpha;
    Beta  = beta;
    
    % Logistic transformation to make the parameters be in the proper range
    Rho    =  -0.99 + (1.98/(1 + exp(-(estpar(1) + estpar(2)*fvec(end)*(fvec(end) <= 0) +...
        estpar(3)*fvec(end)*(fvec(end) > 0) + estpar(4)*rho(end)))));
    Zeta   =  0.1 + (104.9/(1 + exp(-(estpar(5) + estpar(6)*abs(fvec(end))*(fvec(end) <= 0) +...
        estpar(7)*abs(fvec(end))*(fvec(end) > 0) + estpar(8)*(zeta(end))))));

    Alpha  = [Alpha,(((Zeta))/((1 - Rho^2)^(1/2)))];
    Beta   = [Beta,Alpha(end)*Rho];
      
    AlphaMat(:,j)        = Alpha(2:end)';
    BetaMat(:,j)         = Beta(2:end)';
    
    RhoMat(:,j)          = (BetaMat(:,j))./(AlphaMat(:,j));
    ZetaMat(:,j)         = sqrt(fvol(1:end-1).*( (AlphaMat(:,j).^2) - BetaMat(:,j).^2 ));
end

GammaMat = sqrt(AlphaMat.^2 - BetaMat.^2);
Skew     = 3.*BetaMat./(AlphaMat.*sqrt(GammaMat));
Kurt     = 3 + (3*(1+4.*((BetaMat.^2)./(AlphaMat.^2))))./(GammaMat);


if size(AA,1) == size(AA,2)
    FuncVal  = sum(FVal) + T*log(abs(det(AA^(-1))));
else
    FuncVal  = sum(FVal);
end

PF = HtForec^(0.5);
PL = HtLast^(0.5);


if isforec
    VarMat = HtForec;
    P      = PF;
    skew   = Skew(end,:);
    kurt   = Kurt(end,:);
else
    VarMat = HtLast;
    P      = PL;
    skew   = Skew(end-1,:);
    kurt   = Kurt(end-1,:);
end

MM3     = [];
for i=1:size(P,2)
    S=[];
    for j=1:size(P,2)
        for k=1:size(P,2)
            u = 0;
            for t=1:size(P,2) 
                u = u + ((P(i,t))*(P(j,t))*(P(k,t)))*skew(t);
            end
            S(j,k) = u;
        end
    end
    MM3=[MM3 S];
end


MM4=[];
for i=1:size(P,2)
    for j=1:size(P,2)
        S=[];
        for k=1:size(P,2)
            for m=1:size(P,2)
                u = 0;
                for t=1:size(P,2)
                    u = u + (P(i,t))*(P(j,t))*(P(k,t))*(P(m,t))*kurt(t);
                end
                S(k,m)=u;
            end
        end
        MM4=[MM4 S];
    end
end

MMM4=[];
for i=1:size(P,2)
    for j=1:size(P,2)
        S=[];
        for k=1:size(P,2)
            for m=1:size(P,2)
                u = 0;
                for t=1:size(P,2)
                    v = 0;
                    for mm = 1:size(P,2)
                        if mm ~= t
                            u = u + ((P(i,t))*(P(j,t))*(P(k,mm))*(P(m,mm))) + ((P(i,t))*(P(j,mm))*(P(k,t))*(P(m,mm))) +...
                                ((P(i,mm))*(P(j,t))*(P(k,t))*(P(m,mm)));    
                        end
                        v = v + u;
                    end
                    S(k,m)=v;
                end
            end
        end
        MMM4=[MMM4 S];
    end
end

MM4 = MM4 + MMM4;

M2  = AA*VarMat*AA';
M3  = AA*MM3*kron(AA,AA)';
M4  = AA*MM4*kron(kron(AA,AA),AA)';

skkurt = {Skew,Kurt};

end


function [LL2,rho,zeta,alpha,beta,delta,mu] = secstepLL(f,param,paraminit,fvola)

[T,~]   = size(f); 

rho   = [paraminit(1),zeros(1,T-1)];
zeta  = [paraminit(2),zeros(1,T-1)];
alpha = zeros(1,T);
beta  = zeros(1,T);
delta = zeros(1,T);
mu    = zeros(1,T);
ll    = zeros(1,T);
gh    = zeros(1,T);
sig   = zeros(1,T);

for i = 2:T
    
    sig(i) = fvola(i);
    
    rho(i)    =  -0.99 + (1.98/(1 + exp(-(param(1) + param(2)*f(i-1)*(f(i-1) <= 0) + param(3)*f(i-1)*(f(i-1) > 0) + param(4)*rho(i-1)))));
    zeta(i)   =  0.1 + (104.9/(1 + exp(-(param(5) + param(6)*abs(f(i-1))*(f(i-1) <= 0) + param(7)*abs(f(i-1))*(f(i-1) > 0) +...
        param(8)*(zeta(i-1))))));
    alpha(i)  = (((zeta(i)))/((1 - rho(i)^2)^(1/2)));
    beta(i)   = alpha(i)*rho(i);
    delta(i)  = sqrt(sig(i));
    mu(i)     = 0;         
    
    gh(i)     = ((alpha(i)/sqrt(sig(i)))/(pi*(delta(i)/sqrt(sig(i)))))*exp( sqrt((alpha(i)/sqrt(sig(i)))^2 -...
        (beta(i)/sqrt(sig(i)))^2) + (beta(i)/sqrt(sig(i)))*((f(i)-(mu(i)/sqrt(sig(i))))/(delta(i)/sqrt(sig(i)))) )*((sqrt(1+...
        ((f(i)-(mu(i)/sqrt(sig(i))))/(delta(i)/sqrt(sig(i))))^2))^(-1))*besselk(1,(alpha(i)/sqrt(sig(i)))*...
        (sqrt(1+((f(i)-(mu(i)/sqrt(sig(i))))/(delta(i)/sqrt(sig(i))))^2)));
    
    if isnan(gh(i)) || gh(i) == 0
        gh(i) = 0.001;
    end
    
    ll(i)  = log(gh(i));

end

ll  = ll(2:end);

LL2 = -sum(ll);

end

