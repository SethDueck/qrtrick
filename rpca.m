function [ L, S, Y, LRank, U, S_reduced ] = rpca ( X, convergenceFactor, eigShrinkMagnitude, lambda, S0, Y0, LRank0 )
% Decompose data into low-rank part L and sparse S, with error tolerance 
% FALSE: This returns X = L1 + S1 + Y/mu 
% mu controls step size for stepping data into sparse matrix 
% lambda controls maximum eigenvalue (in frobenius norm) 
% Use invalid values to get default mu and lambda 
% From Brunton data book page 124 

doplotting = false; 

% If user-supplied 
if( ~( (1==numel(eigShrinkMagnitude)) && isfinite(eigShrinkMagnitude) && ...
        (1==numel(lambda)) && isfinite(lambda) ) )
    [eigShrinkMagnitude_default, lambda_default] = getDefaultRpcaParams( X ); 
    
    if( ~((1==numel(eigShrinkMagnitude)) && isfinite(eigShrinkMagnitude)) )
        eigShrinkMagnitude = eigShrinkMagnitude_default; 
    end 
    
    if( ~((1==numel(lambda)) && isfinite(lambda)) )
        lambda = lambda_default; 
    end 
end 

thresh = convergenceFactor*norm(X,'fro'); 

L = zeros(size(X)); 
if( all(size(S0)==size(X)) ) 
    S = S0; 
else 
    S = zeros(size(X)); 
end 
if( all(size(Y0)==size(X)) ) 
    Y = Y0; 
else 
    Y = zeros(size(X)); 
end 
if( 1==numel(LRank0) && LRank0>=0 ) 
    LRank = LRank0;
else 
    LRank = 0; 
end 
errork = inf(1); 
count = 0; 

if(doplotting) 
figure(1) 
clf 
hold all 
end 

slist = {}; 
ylist = {}; 
while( (thresh<errork) && (count<1000)); 
    [L, LRank, U, S_reduced] = svt(X-S+eigShrinkMagnitude*Y, eigShrinkMagnitude, max(LRank-1,0)); 
%     [L, LRank] = svt(X-S+(1/mu)*Y, 1/mu, LRank); 
    S = shrink(X-L+eigShrinkMagnitude*Y, lambda*eigShrinkMagnitude); 
    errork = norm(X-L-S,'fro'); 
    slist{count+1} = S; 
    
    if(doplotting) 
    sfac = linspace(0.9,1.1,101);
    serr = nan(size(sfac));
    for sidx = 1:numel(sfac)
        Stemp = shrink(X-L+sfac(sidx)*eigShrinkMagnitude*Y, lambda*eigShrinkMagnitude);
        serr(sidx) = norm(X-L-Stemp,'fro');
    end
    plot(sfac,serr,'.')
    end 
%     figure(1)
%     clf 
%     hold all 
%     temp = (mu*(X-L-S)); 
%     plot(temp(:),Y(:),'.')    
%     title([num2str(errork) '   ' num2str(thresh)] ) 

    Y = Y + (X-L-S)/eigShrinkMagnitude;
    ylist{count+1} = Y; 
    
    count = count + 1; 
end 

nan(0); 

end

