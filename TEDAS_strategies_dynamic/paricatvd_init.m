%{
Function to estimate the initial parameters for ICATVD using the 
GARCH assumption about the independent factors and 
the Normal Inverse Gaussian distribution 
%}
function [FuncVal,F,A,FullPar,InitPar,SigMat] = paricatvd_init(Data,method)

if nargin==2
    if ~ischar(method)
        error('The input "method" must be one of the listed ICA nonlinearity methods.')
    elseif strcmpi(method,'pow3') == 1 || strcmpi(method,'tanh') == 1 || strcmpi(method,'gauss')...
            == 1 || strcmpi(method,'skew') == 1 
        meth = method;       
    else
        error('Unrecognized "method" input.')      
    end
else
    meth = 'pow3';
end

[T,~]         = size(Data);
[FactMat,A,~] = fastica(Data','maxNumIterations',200000,'g',meth);
F             = FactMat';
[~,k]         = size(F);


options  = optimoptions('fmincon','Algorithm','interior-point','TolFun',1e-6,'MaxFunEvals',50000,'MaxIter',20000);

parinit = [1,0.5,1,0.5,0.01,0.8,0.1];
EstPar  = zeros(k,length(parinit)); 
FVal    = zeros(1,k);
AA      = [ 0 0 -1  0  0  0  0;
            0 0 -1 -1  0  0  0;
            0 0 -1  1  0  0  0; 
            0 0  0  0 -1  0  0; 
            0 0  0  0  0 -1  0;
            0 0  0  0  0  0 -1;
            0 0  0  0  1  1  1];
bb      = [0, 0, 0, -1e-10, 0, 0, 0.9999999];

SigMat = zeros(T,k); 

parfor j = 1:k
    facvec = F(:,j);
          
    try
        [estpar,fval] = fmincon(@(param)secdstepLL(facvec,param),parinit,AA,bb,[],[],[],[],[],options);
    catch err
        if any(abs(facvec) > std(facvec)*3)
            ind         = find(abs(facvec) > std(facvec)*3);
            facvec(ind) = mean(facvec);
        end 
        [estpar,fval] = fmincon(@(param)secdstepLL(facvec,param),parinit,AA,bb,[],[],[],[],[],options);
    end
            
    EstPar(j,:)   = estpar;
    FVal(j)       = -fval;
    [~,sig]       = secdstepLL(facvec,estpar);
    sigF          = estpar(5) + estpar(6)*sig(end) + estpar(7)*(facvec(end)^2);
    SigMat(:,j)   = [sig';sigF];
end

RhoMat  = EstPar(:,2)./EstPar(:,1);
ZetaMat = sqrt(( EstPar(:,1)./SigMat(2,:)' ).^2 + ( EstPar(:,2)./SigMat(2,:)' ).^2);
FullPar = EstPar;
InitPar = [RhoMat,ZetaMat];

FuncVal = sum(FVal);

end


function [LL2,sig] = secdstepLL(fac,param)

% Pre-sampling
ifac   = fac(2:end); 
ifac0  = fac(1); 
[T,~]  = size(ifac);

gh     = zeros(1,T);
sig    = zeros(1,T);
delta  = zeros(1,T);
mu     = zeros(1,T);
ll     = zeros(1,T);
rho    = zeros(1,T);

ar0    = param(1);
ar1    = param(2);
alpha  = param(3);
beta   = param(4);
omega  = param(5);
grch   = param(6);
arch   = param(7);


% recursively define the likelihood function with pre-sampling
sig(1) = var(ifac); 
sig(2) = omega + grch*sig(1) + arch*( (ifac(1) - ar0 - ar1*ifac0)^2 );

for i = 3:T  
    
    sig(i)   = omega + grch*sig(i-1) + arch*( (ifac(i-1) - ar0 - ar1*ifac(i-2))^2 );
    
    delta(i) = 1;
    mu(i)    = 0;   
    gh(i)    = ((alpha/sqrt(sig(i)))/(pi*delta(i)))*exp( sqrt((alpha/sqrt(sig(i)))^2 - (beta/sqrt(sig(i)))^2) +...
        (beta/sqrt(sig(i)))*((ifac(i)-mu(i))/delta(i)) )*((sqrt(1+...
        ((ifac(i)-mu(i))/delta(i))^2))^(-1))*besselk(1,(alpha/sqrt(sig(i)))*(sqrt(1+((ifac(i)-mu(i))/delta(i))^2)));           
    ll(i)    = log(gh(i));

end

ll  = ll(3:end);
LL2 = -sum(ll);

end

