function [ out ] = unshrink( X,tau )
% Inverse of shrink, X==shrink(unshrink(X,tau),tau) 

out = sign(X).*(abs(X)+tau); 

end

