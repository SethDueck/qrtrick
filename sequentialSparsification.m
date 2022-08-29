
% Use a "typical" normalized ssalbedo dataset for prototyping 
% ssdata = data; 
if(~exist('ssdata','var')) 
    temp = load('ssdataset.mat'); 
    ssdata = temp.data; 
%     warning('Haven''t acutally loaded this dataset before, probably have to pull data out of a temp var.'); 
    
end 

% Normalize data 
ssdata = ssdata -  repmat(mean(ssdata,  1), [size(ssdata,1),1]); 
ssdata = ssdata ./ repmat(std( ssdata,0,1), [size(ssdata,1),1]); 

% Set options for each step 
srdopts.fidx_iterationConvergence = 111; 
% srdopts.fidx_iterationConvergence = nan(1); 
% errorSequence = 10.^linspace(-2,-5,24); 
errorSequence = 10.^linspace(-3,-3,24); 
convergenceLSError = 1e-3; 
filter_todo = 1:size(ssdata,1); 

% Pre-filter out the worst wavelengths 
toCluster = ssdata(filter_todo,:); 
mtc = mean(toCluster,1); 
dist2_mtc = sum((toCluster-repmat(mtc,[size(toCluster,1) 1])).^2,2); 
[~,leastMeanIdx] = max(dist2_mtc); 
tcml = toCluster-repmat(toCluster(leastMeanIdx,:),[size(toCluster,1) 1]);
dist2_lmi = sum(tcml.^2,2); 
filter_worst = filter_todo(dist2_mtc>findKnee(dist2_mtc)); 
filter_todo(dist2_mtc>findKnee(dist2_mtc)) = []; 
    
numRounds = 24; 
maxLRank = 90; 
subTerm = cell([numRounds 1]); 
Lall = cell([numRounds 1]);  
LRankAll = nan([numRounds 1]); 
Sall = cell([numRounds 1]); 
cindices_clustered = cell([numRounds 1]); 
cindices_fits = cell([numRounds 1]); 
for nidx = 1:numRounds 
    
    if(4000>numel(filter_todo)) 
        break; 
    end 
    
    if(5==nidx) 
        nan(0); 
    end 
    % Recenter data 
    toCluster = ssdata(filter_todo,:); 
    subTerm{nidx} = mean(toCluster,1); 
%     mtc = mean(toCluster,1); 
%     dist2_mtc = sum((toCluster-repmat(mtc,[size(toCluster,1) 1])).^2,2); 
%     [~,leastMeanIdx] = max(dist2_mtc); 
%     tcml = toCluster-repmat(toCluster(leastMeanIdx,:),[size(toCluster,1) 1]);
%     dist2_lmi = sum(tcml.^2,2); 
    
    subTerm{nidx} = toCluster(1,:);
    toCluster = toCluster - repmat(subTerm{nidx},[size(toCluster,1) 1]); 
    
    % Decompose data into L+S, save result 
    [filter_wasClustered, L,LRank, S, timedata] = ...
    sequentialSparsification_step( toCluster, errorSequence, maxLRank, convergenceLSError, srdopts ); 
    Lall{nidx} = L; 
    LRankAll(nidx) = LRank; 
    Sall{nidx} = S; 
    cindices_clustered{nidx} = filter_todo(filter_wasClustered); 
    
    % Do svd on this cluster, include wavelengths that are "as
    % 3-dimensional" as the knee of what was clustered 
    [v_c,d_c] = eig(toCluster(filter_wasClustered,:)'*toCluster(filter_wasClustered,:));
    [d_c,sortfilter] = sort(real(diag(d_c)),'descend'); 
    v_c = v_c(:,sortfilter); 
    p_c = v_c'*toCluster'; 
    non3dPart_mag = sum(p_c(4:end,:).^2); 
    filter_fitsInCluster = non3dPart_mag < findKnee(non3dPart_mag(filter_wasClustered)); 
    cindices_fits{nidx} = filter_todo(filter_fitsInCluster); 
    
%     if(1<nidx)
    sequentialSparsification_cleanManifold; 
%     end
    
    filter_wasClustered_b = false(size(filter_fitsInCluster)); 
    filter_wasClustered_b(filter_wasClustered) = true; 
    
    % Clear clustered wavelengths from todo list 
%     filter_todo(filter_wasClustered) = []; 
    filter_todo(filter_fitsInCluster | filter_wasClustered_b) = []; 
end 

cindices_final = sequentialSparsification_lastClusters( ssdata, filter_todo );
cindices_worst = sequentialSparsification_lastClusters( ssdata, filter_worst ); 

% Break apart any subgroups that were clustered together 
subGroups = cell([numRounds 1]); 
for nidx = 1:numRounds 
cdata = ssdata(cindices_clustered{nidx},:);
cdata = cdata ./ repmat(sqrt(sum(cdata.^2,2)),[1 size(cdata,2)]);
[ut,st,vt] = svd(cdata,'econ');
temp = ut*st;
[~,subgroupScore] = max(abs(temp),[],2);
subGroups{nidx} = subgroupScore; 
end 

Vall = cell([numRounds 1]); 
Pall = cell([numRounds 1]); 
pscores = nan([numRounds size(ssdata,1)]); 
pscores_exp = nan(size(pscores)); 
% Calculate score of each wavelength in each cluster 
for nidx = 1:numel(subTerm) 
    sDataCentered = ssdata - repmat(subTerm{nidx},[size(ssdata,1) 1]); 
    [v,d] = eig(Lall{nidx}'*Lall{nidx}); 
    [~,sortfilter] = sort(diag(d)); 
    v = v(:,sortfilter); 
    Vall{nidx} = v(:,end:-1:(end-LRankAll(nidx)+1));
    p = Vall{nidx}'*sDataCentered'; 
    Pall{nidx} = p; 
    unp = v(:,end:-1:(end-LRankAll(nidx)+1))*p;
    pscores(nidx,:) = sqrt(sum((unp-sDataCentered').^2,1)); 
    pscores_exp(nidx,:) = exp(-5*pscores(nidx,:)/max(pscores(nidx,:)));
end 

numQrPoints = 100; 
qrIndices = nan([numQrPoints 1]); 
qrMag     = nan([numQrPoints 1]); 
qrweights = pscores_exp; 
for nidx = 1:numQrPoints 
    [~,pivotidx] = max(sum(qrweights.^2,1)); 
    qrIndices(nidx) = pivotidx; 
    pivotUnit = qrweights(:,pivotidx); 
    qrMag(nidx) = sqrt(sum(pivotUnit.^2)); 
    pivotUnit = pivotUnit / qrMag(nidx); 
    pUnitWeights = pivotUnit'*qrweights; 
    qrweights = qrweights - pivotUnit*pUnitWeights; 
end 
    

[filter_wasClustered2, L,S, timedata] = ...
    sequentialSparsification_step( ssdata(filter_todo,:), errorSequence, convergenceLSError, srdopts ); 
filter_todo(filter_wasClustered2) = []; 
[filter_wasClustered3, L,S, timedata] = ...
    sequentialSparsification_step( ssdata(filter_todo,:), errorSequence, convergenceLSError, srdopts ); 
filter_todo(filter_wasClustered3) = []; 
[filter_wasClustered4, L,S, timedata] = ...
    sequentialSparsification_step( ssdata(filter_todo,:), errorSequence, convergenceLSError, srdopts ); 
filter_todo(filter_wasClustered4) = []; 


[u_L1,s_L1,v_L1] = svd(L1,'econ'); 
[u_S1,s_S1,v_S1] = svd(S1,'econ'); 
[u_r1,s_r1,v_r1] = svd(ssdata(filter1,:)-L1,'econ'); 

p_L1 = v_L1'*L1'; 
p_S1 = v_S1'*S1'; 
p_r1 = v_r1'*(ssdata(filter1,:)-L1)'; 


