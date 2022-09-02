function [ coeffs ] = qrtrick_invert( basis, target )
% Perform sr3 to get target = basis*wk 

A = basis; 
b = target; 
C = eye(size(A,2)); 
pointWeights = ones(size(b)); 

eta = 1e2; 
lambda = 1e-1; 
w0 = ones([size(A,2) 1]); 
convergenceThreshold = 1e-5; 
maxNumIterations = 1e4; 
[ xk, wk ] = srdsr3_mod( A, C, b, pointWeights, eta, lambda, w0, convergenceThreshold, maxNumIterations ); 

coeffs = xk; 

nan(0); 
    

end

