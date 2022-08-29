function [ out, LRank, U, S_reduced ] = svt( X, tau, maxLRank )
% Singular value threshold operator 
% Takes SVD and throws away small principal components 
% minLRank overrides tau if shrinking to tau would result in
% rank(out)<minLRank 
% From Brunton data book page 124 

% [U,S,V] = svd(X,'econ'); 
[U,S,V] = srd_svd(X,0); 

% Depart from Brunton, allow user to set minimum number of principal
% components to keep -- convergence of rpca suffers if too many pcs are 
% removed at once. 
% slist = diag(S); 
slist = S; 
if( (int8(maxLRank)==maxLRank) && maxLRank>=0 )
    if( 0<maxLRank ) 
%         if( slist(maxLRank) > tau ) 
%             tau = slist(maxLRank+1); 
        if( slist(maxLRank) <= tau ) 
            tau = slist(maxLRank+1); 
        end 
    end 
else 
    warning('Enter non-negative integer for minLRank') 
end 

LRank = sum(slist>tau); 
% S_reduced = shrink(S,tau); 
S_reduced = diag(shrink(S,tau)); 
out = U(:,1:LRank)*(S_reduced(1:LRank,1:LRank)*V(:,1:LRank)'); 
% S(S<tau) = 0; 
% out = U*S*V'; 

end

