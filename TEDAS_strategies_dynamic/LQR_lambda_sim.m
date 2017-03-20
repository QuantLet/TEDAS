function [ final_lambda ] = LQR_lambda_sim( X, tau, numsim, n )

numsim = max(numsim, n);

[ NumRows, NumCols ] = size( X );

U      = rand(NumRows, numsim);
X_norm = zeros(NumCols, 1);
lambda = zeros(numsim, 1);

for j = 1:NumCols
    X_norm(j) = norm(X(:,j));
end    

for k = 1:numsim     
    lambda(k) = max( abs((X'*( (U(:,k)<tau) - tau ) )./X_norm) );        
end

final_lambda = quantile(lambda, 0.1); 

end