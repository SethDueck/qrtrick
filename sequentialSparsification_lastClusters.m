function [ cindices_final] = sequentialSparsification_lastClusters( ssdata, filter_todo )


% Now use google method to cluster remaining wavelengths 
pmatrix = (exp(-squareform(pdist(ssdata(filter_todo,:)))));
for yidx = 1:size(pmatrix,2)
    pmatrix(:,yidx) = pmatrix(:,yidx)/sqrt(sum(pmatrix(:,yidx).^2));
end
[vp,dp] = eigs(pmatrix,10);
p = vp'*pmatrix;
[va,da] = eig(abs(p)*abs(p)');
[da,sortfilter] = sort(real(diag(da)),'descend'); 
va = va(:,sortfilter); 
numFinalGroups = sum(da > (sum(da)/10)); 
[~,dissimilarPDims] = max(abs(va(:,1:numFinalGroups)),[],2); 
finalGroupScores = nan([numel(filter_todo) numFinalGroups]); 
for nidx = 1:numFinalGroups 
    finalGroupScores(:,nidx) = sum(p(dissimilarPDims==nidx,:).^2,1); 
end 

% [~,finalGroupIdx] = max(finalGroupScores,[],2); 
[~,finalGroupIdx] = max(dp*abs(p),[],1);

cindices_final = cell([max(finalGroupIdx) 1]); 
for nidx = 1:numel(cindices_final) 
    cindices_final{nidx} = filter_todo(finalGroupIdx==nidx); 
end 



end 

