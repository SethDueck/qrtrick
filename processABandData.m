
% Load results from a-band-test data -- it's #test run on the last
% wavelength 
load('abandTest.mat'); 

% Normalize elements of Z to their value at lowest scale factor 
vectorShapes = nan([size(mmFactors,3) size(mmFactors,1)*size(mmFactors,2)]); 
mmVector = nan(size(vectorShapes)); 
daVector = nan(size(vectorShapes)); 
runningIdx = 0; 
for didx = 1:size(mmFactors,1) 
    sumationTerms = (squeeze(mmFactors(didx,:,:).*daFactors(didx,:,:)))'; % matrix is nDims by nDims 
    summationShape = sumationTerms./repmat(sumationTerms(1,:),[size(sumationTerms,1) 1]); 
    vectorShapes(:,runningIdx+(1:size(mmFactors,1))) = summationShape; 
    mmVector(:,runningIdx+(1:size(mmFactors,1))) = squeeze(mmFactors(didx,:,:))';
    daVector(:,runningIdx+(1:size(mmFactors,1))) = squeeze(daFactors(didx,:,:))';
    runningIdx = runningIdx + size(mmFactors,1); 
end 
vectorShapes_scaled = vectorShapes .* repmat(scalefactors',[1 size(vectorShapes,2)]); 
vectorShapes_norm = vectorShapes_scaled ./ repmat(sqrt(sum(vectorShapes_scaled.^2,1)),[size(vectorShapes_scaled,1) 1]); 
distanceFromOthers = (1-(vectorShapes_norm'*vectorShapes_norm));
sumDistanceFromOthers = sqrt(sum(distanceFromOthers)); 
filter_isBestBehaved = sumDistanceFromOthers < 1.21*min(sumDistanceFromOthers); 
filter_isBestBehaved = abs(sumDistanceFromOthers-min(sumDistanceFromOthers)) < 0.2*std(sumDistanceFromOthers); 
bestBehavedShape = mean(vectorShapes_norm(:,filter_isBestBehaved),2); 
oddlyBehavedShape = vectorShapes_norm - repmat(bestBehavedShape,[1 size(vectorShapes,2)]); 
[v,d] = eig(oddlyBehavedShape*oddlyBehavedShape');
p2 = oddlyBehavedShape'*v; 

temp = vectorShapes;
temp = temp -  repmat(mean(temp,2),  [1 size(temp,2)]);
temp = temp ./ repmat(std( temp,0,2),[1 size(temp,2)]);
temp(isnan(temp)) = 0;
[v,d] = eig(temp*temp'); 
p2 = temp'*v;

[v_mm,d_mm] = eig(mmVector*mmVector'); 
[v_da,d_da] = eig(daVector*daVector'); 

% It looks like the "best" (regular) factor profiles all have the same
% "shape" of #mmFactor. 
mm_norm1 = mmVector./repmat(mmVector(1,:),[size(mmVector,1) 1]); 
mm_norm1 = mm_norm1./repmat(sqrt(sum(mm_norm1.^2,1)),[size(mm_norm1,1) 1]); 
% mm_norm1 = mmVector; 
da_norm1 = daVector./repmat(daVector(1,:),[size(daVector,1) 1]); 
doCompareMMFactors = true; 
if( doCompareMMFactors )
figure(19)
clf
subplot(2,2,1)
surf(mm_norm1(:,filter_isBestBehaved)); 
view([0 0 1])
title('mmVector')
ylabel('best')
subplot(2,2,2)
surf(da_norm1(:,filter_isBestBehaved)); 
view([0 0 1])
title('daVector')
subplot(2,2,3)
surf(mm_norm1(:,~filter_isBestBehaved)); 
view([0 0 1])
ylabel('~best')
subplot(2,2,4)
surf(da_norm1(:,~filter_isBestBehaved)); 
view([0 0 1])
end 
% 
%
% Now isolate the "anomaly" in the mmFactor and try to predict the anomaly
% from the data. 
%
bestBehavedMM = mean(mm_norm1(:,filter_isBestBehaved), 2); 
bestBehavedMM = bestBehavedMM / sqrt(sum(bestBehavedMM.^2)); % Normalize 
% badMMs = mmVector(:,~filter_isBestBehaved)./repmat(mmVector(1,~filter_isBestBehaved),[size(daVector,1) 1]); 
goodComponentMagnitude = bestBehavedMM'*mm_norm1; 
badMMs_badComponent = mm_norm1 - bestBehavedMM*goodComponentMagnitude;
[v_badMMComponents,d_badMMComponents] = eig(badMMs_badComponent*badMMs_badComponent'); 
p_badMMComponents = badMMs_badComponent'*v_badMMComponents; 
% p_badMMComponents(:,end-0) can be taken as the "score" of how the MM
% deviates from the mean of the "good" MMs 

