function [ referencePixelIndices ] = srd_findMostInformativePixels( pdistMatrix, numToFind )

data_mean = mean(pdistMatrix,1); 
data_norm = pdistMatrix -  repmat(data_mean,[size(pdistMatrix,1) 1]); 
% rid_mean = mean(rid_norm,2); 
% rid_norm = rid_norm -  repmat(rid_mean,[1 size(riDistances,1)]); 
% Columns have zero mean -- are "extraordinary" distance from this point to
% point rowidx 
[v,~] = eig(data_norm' * data_norm  ); 
% Decide how many eigen-distance-functions to keep 

[~,filter] = sort(sum(abs(v),1));
v_toSearch = v(:,filter(end-numToFind+1:end)); 
% vrid_toSearch = vrid; 
% Find most informative pixels to sample 

referencePixelIndices = nan([numToFind 1]); 
thisBestPixel = nan(1);         
for vidx = 1:numToFind 
    if( 1 < vidx ) 
        % Orthogonalize eigenfunctions wrt pixels in sample list 
        lastBestPixel = v_toSearch(thisBestPixel,:) / sqrt(sum(v_toSearch(thisBestPixel,:).^2)); 
        v_toSearch = v_toSearch - ...
            (v_toSearch*lastBestPixel')*lastBestPixel; 
    end 
    [~,thisBestPixel] = max(sum(v_toSearch.^2,2),[],1); 
    referencePixelIndices(vidx) = thisBestPixel; 
end 



end

