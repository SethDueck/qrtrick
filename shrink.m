function [ out ] = shrink( X,tau )
% From Brunton book page 124 

out = sign(X).*max(abs(X)-tau,0); 

end

