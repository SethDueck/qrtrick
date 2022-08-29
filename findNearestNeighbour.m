function [ dnn ] = findNearestNeighbour( A, B )
% A, B are nPoints-by-nDims 

ndims = size(A,2); 

isConverged = false; 
Aidx_old = nan(1); 
Bidx_old = nan(1); 
Ann = mean(A,1); 
Bnn = mean(B,1); 
dnn = sum((Ann-Bnn).^2); 
flagAB = true(1); 

while ~isConverged 

if( flagAB )
    % Move A point closer 
    temp = A; 
    for didx = 1:ndims 
        temp(:,didx) = temp(:,didx)-Bnn(didx); 
    end 
    ds = sum(temp.*temp,2); 
    [dnn,Aidx] = min(ds); 
    if( Aidx_old==Aidx ) 
        isConverged = true; 
    else 
        Aidx_old = Aidx; 
        Ann = A(Aidx,:); 
    end 
else 
    % Move B point closer 
    temp = B; 
    for didx = 1:ndims 
        temp(:,didx) = temp(:,didx)-Ann(didx); 
    end 
    ds = sum(temp.*temp,2); 
    [dnn,Bidx] = min(ds); 
    if( Bidx_old==Bidx ) 
        isConverged = true; 
    else 
        Bidx_old = Bidx; 
        Bnn = B(Bidx,:); 
    end 
end 

flagAB = ~flagAB; 

end 

dnn = sqrt(dnn); 

end

