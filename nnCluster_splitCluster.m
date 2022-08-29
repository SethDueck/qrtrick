function [ subclusterIdx, methodsToTry ] = nnCluster_splitCluster( p, thisCluster, methodsToTry ) 

keepTryingSplit = true; 
subclusterIdx = ones([size(p,2) 1]); 
nextSubIdx = 2; 

if( methodsToTry.gaps && keepTryingSplit ) 
% Look for gaps in the data 
for didx = 1:size(p,1) 
    
        [un,ux] = hist(p(didx,:)); 
        [~,midx] = min(un); 
        utSplitPoint = ux(midx); 
%         uts = sort(ut(:,uidx)); 
%         [~,midx] = max(diff(uts)); 
%         utSplitPoint = 0.5*(uts(midx)+uts(midx+1)); 
        std_whole = std(p(didx,:)); 
        filter_gt = p(didx,:) >= utSplitPoint; 
        std_pos = std(p(didx, filter_gt)); 
        std_neg = std(p(didx,~filter_gt)); 
        if( std_whole > 2*sqrt(std_pos*std_pos+std_neg*std_neg) ) 
            % Do split; send smaller group to a new cluster 
            keepTryingSplit = false; % Split has been achieved 
            if(sum(filter_gt+0) > 0.5*size(p,2)) 
               subclusterIdx(~filter_gt) = nextSubIdx; 
            else 
               subclusterIdx( filter_gt) = nextSubIdx;  
            end
            
            % Only split one group this this -- usually if this split is
            % triggered it's because there's a single "outlier cluster"
            % that could trigger the break in several dimensions. 
            break; 
        end 
        
end 

if( keepTryingSplit ) 
    % This method has been exhausted, don't try it again  
    methodsToTry.gaps = false; 
end 

end 


if( methodsToTry.tails && keepTryingSplit ) 
    
% Look for bad tails in the data 
findKneeValues = nan([size(p,2) (size(p,1)+(size(p,1)*(size(p,1)-1))/2)]); 
kvidx = 1; 
for didx1 = 1:size(p,1) 
    findKneeValues(:,kvidx) = p(didx1,:); 
    kvidx = kvidx+1; 
    for didx2 = (didx1+1):size(p,1) 
        findKneeValues(:,kvidx) = p(didx1,:)+p(didx2,:); 
        kvidx = kvidx+1; 
    end 
end 

tailRatios_n = nan([size(findKneeValues,2) 1]); 
tailRatios_p = nan(size(tailRatios_n)); 
for kvidx = 1:size(findKneeValues,2) 
    [~,~,badnessStat] = findKnee_subSlope(findKneeValues(:,kvidx)); 
    tailRatios_n(kvidx) = badnessStat(1)/badnessStat(2); 
    tailRatios_p(kvidx) = badnessStat(3)/badnessStat(2); 
end 
% Filter and move "weird tails" to end of clabels list 
tailRatios_n(isnan(tailRatios_n)) = 0; 
tailRatios_p(isnan(tailRatios_p)) = 0; 
[~,findKneeIdx] = max(tailRatios_n+tailRatios_p); 
[kneeval_n,kneeval_p,~] = findKnee_subSlope(findKneeValues(:,findKneeIdx)); 
if(isnan(kneeval_n)) 
    kneeval_n = -inf(1); 
end 
if(isnan(kneeval_p)) 
    kneeval_p = inf(1); 
end 
filter_knees = (findKneeValues(:,findKneeIdx)>kneeval_n) & (findKneeValues(:,findKneeIdx)<kneeval_p); 
subclusterIdx(~filter_knees) = nextSubIdx; 

% Try this method once and then move on 
methodsToTry.tails = false; 
keepTryingSplit = false; 
end 

if( methodsToTry.unfold && keepTryingSplit ) 
    [subgroupIndicator, nongroupIndicator] = nnCluster_splitParallelPlanes( thisCluster ); 

    % Assume subgroupIndicator splits on +/-, but throw away wavelengths
    % that are not well-differentiated 
    filter_isDifferentiated = abs(subgroupIndicator) > nongroupIndicator; 
    if( sum(subgroupIndicator(filter_isDifferentiated)< 0) > ... 
        sum(subgroupIndicator(filter_isDifferentiated)>=0) ) 
        subclusterIdx( subgroupIndicator< 0 ) = nextSubIdx; 
    else 
        subclusterIdx( subgroupIndicator>=0 ) = nextSubIdx; 
    end 
    nextSubIdx = nextSubIdx + 1; 
    subclusterIdx( ~filter_isDifferentiated ) = nextSubIdx; 
    
    % Use this method once, then allow other methods to be used again on
    % the resulting subclusters 
    methodsToTry.unfold = false;
    methodsToTry.gaps = true; 
    methodsToTry.tails = true; 
    keepTryingSplit = false; 
end 


% [ kneeVal_n, kneeVal_p, fitBadnessMeasure ] = findKnee_subSlope( p(didx,:) );  
%     if(fitBadnessMeasure(1)/fitBadnessMeasure(2) > 15) 
%         % Split on kneeVal_n 
%         subclusterIdx(p(didx,:)<kneeVal_n) = nextSubIdx; % Move "bad" cluster to the end 
%         nextSubIdx = nextSubIdx + 1; 
%     end 
%     if(fitBadnessMeasure(3)/fitBadnessMeasure(2) > 15) 
%         % Split on kneeVal_p 
%         subclusterIdx(p(didx,:)>kneeVal_p) = nextSubIdx; % Move "bad" cluster to the end 
%         nextSubIdx = nextSubIdx + 1;         
%     end 
%  
% end 



end 



