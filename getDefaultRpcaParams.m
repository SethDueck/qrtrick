function [eigShrinkMagnitude, lambda_default] = getDefaultRpcaParams( X ) 

[n1,n2] = size(X); 
eigShrinkMagnitude = 4*sum(abs(X(:))) / (1*n1*n2); 
lambda_default = 1/sqrt(max(n1,n2)); 

end

