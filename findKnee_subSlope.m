function [ kneeVal_n, kneeVal_p, fitBadnessMeasure ] = findKnee_subSlope( data )

% Find knee on high end of data 
data_sorted = sort(data(:)); 
data_xFrom1  = linspace(1,0,numel(data_sorted))'; 
data_y  = (data_sorted-data_sorted(1)) / (data_sorted(end)-data_sorted(1)); 
data_y = data_y - median(data_y); 
data_r2 = data_xFrom1.^2 + data_y.^2; 
[~,data_idxOfSortedKnee_minr2] = min(data_r2);  
[~,data_idxOfSortedKnee_maxRatio] = max(data_y./data_r2);  
idxOfSortedKnee_p = min(data_idxOfSortedKnee_minr2,data_idxOfSortedKnee_maxRatio);

% Find knee on low end of data 
data_y_flip = -flipud(data_y);
data_r2 = data_xFrom1.^2 + data_y_flip.^2; 
[~,data_idxOfSortedKnee_minr2] = min(data_r2);  
[~,data_idxOfSortedKnee_maxRatio] = max(data_y_flip./data_r2);  
idxOfSortedKnee_n = numel(data) - min(data_idxOfSortedKnee_minr2,data_idxOfSortedKnee_maxRatio) + 1;

% Fit line to data between knees 
% p = polyfit(data_xFrom1(idxOfSortedKnee_n:idxOfSortedKnee_p), data_y(idxOfSortedKnee_n:idxOfSortedKnee_p), 1); 
pointWeights = ones([(idxOfSortedKnee_p-idxOfSortedKnee_n+1) 1]); 
A = [data_xFrom1(idxOfSortedKnee_n:idxOfSortedKnee_p),ones(size(pointWeights))]; 
b = data_y(idxOfSortedKnee_n:idxOfSortedKnee_p); 
C = eye(size(A,2)); 
eta = 1e0;
lambda = 1e-1; 
w0 = ones([size(A,2) 1]); 
convergenceThreshold = 1e-5; 
maxNumIterations = 1e4; 
[ xk, wk ] = srdsr3_mod( A, C, b, pointWeights, eta, lambda, w0, convergenceThreshold, maxNumIterations ); 

data_y_shifted = data_y - polyval(xk',data_xFrom1); 

[~,midx_n_rough] = min(data_y_shifted.^2+flipud(data_xFrom1).^2); 
[~,midx_p_rough] = min(data_y_shifted.^2+(data_xFrom1).^2); 
fit_mean = mean(data_y_shifted(midx_n_rough:midx_p_rough)); 
fit_std  = std( data_y_shifted(midx_n_rough:midx_p_rough)); 
 
midx_n = find(data_y_shifted>(fit_mean-fit_std),1,'first');
midx_p = find(data_y_shifted<(fit_mean+fit_std),1,'last'); 
kneeVal_n = data_sorted(midx_n); 
kneeVal_p = data_sorted(midx_p); 

% Score confidence in knee positions 
m1 = nan(3,1); 
area_total = sum(abs(data_y_shifted)); 
m1(1) = sum(abs(data_y_shifted(1:midx_n-1))); 
m1(2) = sum(abs(data_y_shifted(midx_n:midx_p))); 
m1(3) = sum(abs(data_y_shifted((midx_p+1):end))); 

m2 = nan(3,1); 
area_total = sum(abs(data_y_shifted)); 
m2(1) = sum(data_y_shifted(1:midx_n-1).^2); 
m2(2) = sum(data_y_shifted(midx_n:midx_p).^2); 
m2(3) = sum(data_y_shifted((midx_p+1):end).^2); 

fitBadnessMeasure = m2 ./ m1; 
nan(0); 

end

