function [ kneeVal, data_idxOfSortedKnee ] = findKnee( data )
% Draw sorted data on x=[0,1], y=[0,1]. Return data value of point on curve 
% closest to (x,y) = (1,0). 

data_sorted = sort(data(:)); 
data_xFrom1  = linspace(1,0,numel(data_sorted))'; 
data_y  = (data_sorted-data_sorted(1)) / (data_sorted(end)-data_sorted(1)); 
data_y = data_y - median(data_y); 
data_r2 = data_xFrom1.^2 + data_y.^2; 
[~,data_idxOfSortedKnee_minr2] = min(data_r2);  
[~,data_idxOfSortedKnee_maxRatio] = max(data_y./data_r2);  
data_idxOfSortedKnee = max(data_idxOfSortedKnee_minr2,data_idxOfSortedKnee_maxRatio); 
% [~,data_idxOfSortedKnee] = min(data_r2); 
kneeVal = data_sorted(data_idxOfSortedKnee); 
    
end

