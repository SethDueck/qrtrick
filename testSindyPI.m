

testDataSet = rotatedPThisCentroid; 

testDataSet = testDataSet -  repmat(mean(testDataSet,2),    [1 size(testDataSet,2)]); 
testDataSet = testDataSet ./ repmat( std(testDataSet,[],2), [1 size(testDataSet,2)]); 
convergenceThresholdFactor = 1e-5; 
[L,S] = rpca(testDataSet,convergenceThresholdFactor);
sparseMag = sqrt(sum(S.^2,1)); 
sparseScore = exp(-5 * sparseMag / max(sparseMag));
L_scored = L.*repmat(sparseScore,[size(L,1) 1]);

% filter_wellBehaved = std(testDataSet,[],1) < 4*median(std(testDataSet,[],1)); 
% testDataSet = testDataSet(:,filter_wellBehaved); 
% testDataSet = testDataSet -  repmat(mean(testDataSet,2),    [1 size(testDataSet,2)]); 
% testDataSet = testDataSet ./ repmat( std(testDataSet,[],2), [1 size(testDataSet,2)]); 


% Clear out outliers 
% testDataSet = repmat(exp(-prod(abs(testDataSet),1)),[size(testDataSet,1) 1]);

% Create polynomial for fitting surface  
powersForAllDimensions = 0:2; 
numel(powersForAllDimensions);  
% Produce grid same as what would be provided by ndgrid 
tempx = ones(ones([1 size(testDataSet,1)])*numel(powersForAllDimensions)); 
for padidx = 1:numel(powersForAllDimensions) 
    tempx(padidx,:) = powersForAllDimensions(padidx); 
end 
% Permute grid dimensions to get each of the different outputs from ndgrid 
polyPowers = nan([numel(tempx),size(testDataSet,1)]); 
tempx_dimOrders_orig = 1:numel(size(tempx)); 
tempx_dimOrders = tempx_dimOrders_orig; 
for padidx = 1:size(testDataSet,1) 
    tempx_dimOrders(padidx) = 1; 
    tempdata = permute(tempx,tempx_dimOrders); 
    polyPowers(:,padidx) = tempdata(:); 
    tempx_dimOrders(padidx) = tempx_dimOrders_orig(padidx)+1; 
end 

%  
% [tempx,tempy,tempz,tempq] = ndgrid(0:2); 
% polyPowers = [tempx(:),tempy(:),tempz(:)]; 
% 
% polyPowers = zeros([6 size(testDataSet,1)]); 
% polyPowers(1 ,     :) = 0; 
% polyPowers(2 , end-0) = 1; 
% polyPowers(3 , end-1) = 1; 
% polyPowers(4 , end-0) = 2; 
% polyPowers(5 , end-1) = 2; 
% polyPowers(6 , end-0) = 1; 
% polyPowers(6 , end-1) = 1; 
% polyPowers(7 , end-0) = 3; 
% polyPowers(8 , end-1) = 3; 
% polyPowers(9 , end-0) = 2; 
% polyPowers(9 , end-1) = 1; 
% polyPowers(10, end-0) = 1; 
% polyPowers(10, end-1) = 2; 
% polyPowers(11, end-0) = 4; 
% polyPowers(12, end-1) = 4; 
% polyPowers(13, end-0) = 3; 
% polyPowers(13, end-1) = 1; 
% polyPowers(14, end-0) = 1; 
% polyPowers(14, end-1) = 3; 
% polyPowers(15, end-0) = 2; 
% polyPowers(15, end-1) = 2; 
srdPolyTerms = srdPolyval(testDataSet, polyPowers); 
srdPolyTerms(:,2:end) = srdPolyTerms(:,2:end) -  repmat(mean(srdPolyTerms(:,2:end),1),[size(srdPolyTerms,1) 1]); 
srdPolyTerms(:,2:end) = srdPolyTerms(:,2:end) ./ repmat(std(srdPolyTerms(:,2:end),[],1),[size(srdPolyTerms,1) 1]); 
% srdPolyTerms(:,2:end) = srdPolyTerms(:,2:end) -  repmat(min(srdPolyTerms(:,2:end),[],1),[size(srdPolyTerms,1) 1]); 
% srdPolyTerms(:,2:end) = srdPolyTerms(:,2:end) ./ repmat(max(srdPolyTerms(:,2:end),[],1),[size(srdPolyTerms,1) 1]); 
% srdPolyTerms(:,2:end) = 2*(srdPolyTerms(:,2:end)-0.5); 
LPolyTerms = srdPolyval(L_scored, polyPowers);
LPolyTerms(:,2:end) = LPolyTerms(:,2:end) -  repmat(mean(LPolyTerms(:,2:end),1),[size(LPolyTerms,1) 1]);
LPolyTerms(:,2:end) = LPolyTerms(:,2:end) ./ repmat(std(LPolyTerms(:,2:end),[],1),[size(LPolyTerms,1) 1]);
[vL,dL] = eig(LPolyTerms'*LPolyTerms); 

% Try sindy-pi on every column of polyterms 
allxk = nan([1 1]*size(polyPowers,1)); 
allResiduals = nan(size(srdPolyTerms)); 
for targetIdx = 1:size(srdPolyTerms,2) 
    
% srdtarget = srdPolyTerms(:,targetIdx)';
% srdbasis  = srdPolyTerms; 
srdtarget = LPolyTerms(:,targetIdx)';
srdbasis  = LPolyTerms; 

maxNumElementsInL1 = 800; 
if( maxNumElementsInL1 < size(srdbasis,1) ) 
    % Use up to maxNumElementsInL1 -- take random sample but maintain order
    filter_selectRandomUnique = sort(randperm(size(srdbasis,1), maxNumElementsInL1)); 
    % Use maxNumElementsInL1 elements with smallest sparse component 
    [~,sparseSortIdx] = sort(sparseScore); 
    filter_leastSparseElements = sparseSortIdx(end-maxNumElementsInL1+1:end); 
    % Use sparse score as a probability distribution 
    filter_sparsePDist = floor(interp1(cumsum(sparseScore)/sum(sparseScore),1:numel(sparseScore),rand([maxNumElementsInL1 1]))); 
    % Choose which filter to use 
%     filter_chooseBasisElements = filter_selectRandomUnique; 
%     filter_chooseBasisElements = filter_leastSparseElements; 
    filter_chooseBasisElements = filter_sparsePDist; 
    
    x      = srdbasis(filter_chooseBasisElements,:); 
    srdbeq = srdtarget(filter_chooseBasisElements); 
    f_residualWeights = fitWeightFn(filter_chooseBasisElements); 
else 
    % Use whole dataset 
    x      = srdbasis; 
    srdbeq = srdtarget; 
    f_residualWeights = fitWeightFn; 
end

% Throw away trivial solutions in sindy-pi 
x(:,targetIdx) = 0; 
srdbasis_zeroed = srdbasis; 
srdbasis_zeroed(:,targetIdx) = 0; 

% Set up tuning parameters for sindy 
C = eye(size(x,2)); 
eta = 0.02; 
lambda = 1e-2; 
% lambda = 1e-11; 
% lambda = 1e-1; 
w0 = ones([size(x,2) 1]); 

% Run L1 optimization and evaluate error 
[xk, wk] = srdsr3( x, C, srdbeq', eta, lambda, w0, 1e-2, 1e3 ); 
xk_zeroed = xk; % Take "better" values from xk, but use same sparsity as wk 
xk_zeroed(0==wk) = 0; 
bRetrieved_fitData = (x*xk)';
residuals_fitData = abs(srdbeq-bRetrieved_fitData); 
bRetrieved = (srdbasis_zeroed*xk)'; 
residuals = abs(srdtarget-bRetrieved); 

% Save retrieved weights and residuals for each fit 
allxk(:,targetIdx) = xk'; 
allResiduals(:,targetIdx) = residuals'; 

end 

% Normalize the implicit equation parameters 
allxk_normalized = allxk ./ repmat(sqrt(sum(allxk.^2,1)),[size(allxk,1) 1]); 
