function [ covMatrix, mostLinearDep] = qrtrick_makeDataCovMatrix( data, covMatrix )
% Replace nans in #covMatrix with covariance measure of #data 

for bidx1 = 1:size(data,2) 
    for bidx2 = (bidx1+1):size(data,2)
        if(isnan(covMatrix(bidx1,bidx2)) ) 
            covMatrix(bidx1,bidx2) = det(data(:,[bidx1 bidx2])'*data(:,[bidx1 bidx2]));             
            covMatrix(bidx2,bidx1) = covMatrix(bidx1,bidx2); 
        end             
    end 
end             
bmax = max(max(abs(covMatrix))); 
mostLinearDep = nan([1 1]*size(covMatrix,2)); % Each row is sorted, best candidate for xidx is at yidx=1 
% linearDepMeasure = nan([1 1]*size(dCovMatrix,2)); 
for bidx = 1:size(covMatrix,2) 
    covMatrix(bidx,bidx) = 1.1*bmax; % Set diagonal to large value so we can find minimum of non-zero values 
%     [linearDepMeasure(bidx,:),mld] = sort(dCovMatrix(bidx,:)); 
    [~,mld] = sort(covMatrix(bidx,:)); 
    mostLinearDep(bidx,:) = mld;
end 

end

