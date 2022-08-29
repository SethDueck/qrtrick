function [ subgroupIndicator, nongroupIndicator ] = nnCluster_splitParallelPlanes( thisCluster )

% There are modes of being distant, and when they are uncorrelated with
% pc1, there is a symmetry about the pc1 axis. 

% Center cluster, project onto pcs 
centeredCluster = thisCluster; 
centeredCluster = centeredCluster - repmat(mean(centeredCluster,1),[size(centeredCluster,1) 1]); 
[vc,dc] = eig(centeredCluster'*centeredCluster); 
[dc,sortfilter] = sort(real(diag(dc)),'descend'); 
sortfilter((dc/sum(dc))<1e-3) = []; % Clear out unimportant dimensions 
vc = vc(:,sortfilter); 
pc = vc'*centeredCluster'; 

% Normalize axes 
pcn = pc; 
for didx = 1:size(pc,1) 
    pcn(didx,:) = pcn(didx,:) / std(pcn(didx,:)); 
end 
mdl = KDTreeSearcher(pcn(1:3,:)','BucketSize',10);
[nnIdx,~] = knnsearch(mdl,pcn(1:3,:)','k',20); 

needsGroup = true([size(pcn,2) 1]); 
filter = needsGroup; 
groupIdx = nan(size(needsGroup)); 
nextGroupIdx = 0; 
while sum(needsGroup) > 0.05*numel(needsGroup) || sum(filter)>0.5*size(nnIdx,2) 
    nextGroupIdx = nextGroupIdx + 1; 
    [numNeighboursNeedGroup,needIdx] = max(needsGroup.*sum(needsGroup(nnIdx),2));
%     needIdx = find(needsGroup,1,'first'); 
    filter = needsGroup & any(nnIdx==needIdx,2); 
    needsGroup(filter) = false; 
    groupIdx(filter) = nextGroupIdx; 
end  

% Find centroids for each group 
centroids = nan([size(pcn,1) nextGroupIdx]); 
for gidx = 1:nextGroupIdx 
    centroids(:,gidx) = mean(pcn(:,gidx==groupIdx),2); 
end 

% Get triangulation
k = delaunay(centroids(1:3,:)');
t = turnTriIntoSurface(k); 
t_edges = unique([t(:,[1 2]); t(:,[1 3]); t(:,[2 3])],'rows'); 

% Create connectivity graph for each group 
% connectivityGraph = zeros(nextGroupIdx*[1 1]); 
% for gidx1 = 1:nextGroupIdx 
%     for gidx2 = (gidx1+1):nextGroupIdx 
%         dnn = findNearestNeighbour( pcn(1:3,gidx1==groupIdx)', pcn(1:3,gidx2==groupIdx)' ); 
%         connectivityGraph(gidx1,gidx2) = dnn; 
%         connectivityGraph(gidx2,gidx1) = dnn; 
%     end 
% end 
connectivityGraph = inf(nextGroupIdx*[1 1]); 
for eidx = 1:size(t_edges,1) 
    dnn = findNearestNeighbour( pcn(1:3,t_edges(eidx,1)==groupIdx)', pcn(1:3,t_edges(eidx,2)==groupIdx)' ); 
    connectivityGraph(t_edges(eidx,1),t_edges(eidx,2)) = dnn; 
    connectivityGraph(t_edges(eidx,2),t_edges(eidx,1)) = dnn; 
end 
reducedConnectivityGraph = reduceConnectivityGraph( connectivityGraph ); 
% [v,d] = eig(reducedConnectivityGraph'*reducedConnectivityGraph); 
% [d,sortfilter] = sort(real(diag(d)),'descend');
% v = v(:,sortfilter);
% p = v'*reducedConnectivityGraph';

% Columns of reduced connectivity graph tends to be linear, but has high 
% values for large distances -- invert this 
nearnessGraph = reducedConnectivityGraph; 
for gidx = 1:size(nearnessGraph,2) 
    nearnessGraph(:,gidx) = max(nearnessGraph(:,gidx)) - nearnessGraph(:,gidx); 
    nearnessGraph(gidx,gidx) = 0; 
    nearnessGraph(:,gidx) = nearnessGraph(:,gidx) / sum(nearnessGraph(:,gidx)); 
end 
dualNearnessGraph = min(nearnessGraph,nearnessGraph'); 
[vdn,ddn] = eig(dualNearnessGraph'*dualNearnessGraph); 
[ddn,sortfilter] = sort(real(diag(ddn)),'descend'); 
vdn = vdn(:,sortfilter); 
pdn = vdn'*dualNearnessGraph'; 

% Find fold differentiation dim in projection of connectivity graph. 
% Usually first pc contains background positive signal, fold
% differentiation is contained in second pc, which is two-signed 
% scatter3(pc(3,~isnan(groupIdx)),pc(2,~needsGroup),pc(1,~needsGroup),5,p(foldDiffDim,groupIdx(~needsGroup)))
isTwoSigned = any(pdn>0,2) & any(pdn<0,2); 
foldDiffDim = find(isTwoSigned,1,'first'); 

% Take fold indicator dimension and magnitude of "anything else" dimension 
subgroupIndicatorVals = pdn(foldDiffDim,:); 
nongroupIndicatorVals = sqrt(sum(pdn(foldDiffDim+1:end,:).^2)); 

subgroupIndicator = nan(size(groupIdx)); 
nongroupIndicator = nan(size(groupIdx)); 
for pidx = 1:numel(subgroupIndicator) 
%     if( isnan(groupIdx(pidx)) ) 
        neighbourGroupIdxs = groupIdx(nnIdx(pidx,:)); 
        neighbourIndicators = subgroupIndicatorVals( neighbourGroupIdxs(~isnan(neighbourGroupIdxs))); 
        subgroupIndicator(pidx) = mean(neighbourIndicators); 
        neighbourIndicators = nongroupIndicatorVals( neighbourGroupIdxs(~isnan(neighbourGroupIdxs))); 
        nongroupIndicator(pidx) = mean(neighbourIndicators); 
%     else 
%         subgroupIndicator(pidx) = p(foldDiffDim,groupIdx(pidx)); 
%     end     
end 

nan(0); 


% 
% % Find modes of being distant 
% dmat = squareform(pdist(pc')); 
% [vd,dd] = eig(dmat'*dmat);
% [dd,sortfilter] = sort(real(diag(dd)),'descend'); 
% sortfilter((dd/sum(dd))<1e-3) = []; % Clear out unimportant dimensions 
% vd = vd(:,sortfilter); 
% pd = vd'*dmat'; 
% 
% % Assume fold is somewhat along the axis of first pc, otherwise it could
% % just be fit to a manifold. This is same as saying pc1 is not single-
% % valued in pc2, pc3, ... 
% dmat1 = squareform(pdist(pc(1,:)')); 
% corrWithPc1 = mean(dmat1*pd')./dd(1:5)'; 
% [~,leastCorrIdx] = min(abs(corrWithPc1)); 
% 
% % Points that are "across each other" in the fold will have distance in
% % dmat but not in dmat_lc. 
% dmat_lc = squareform(pdist(pd(leastCorrIdx,:)')); 
% 
% % Normalize the distance matrices 
% dmat_norm    = dmat; 
% dmat_lc_norm = dmat_lc; 
% for didx = 1:size(dmat,2) 
%     dmat_norm(:,didx) = dmat_norm(:,didx) / sum(dmat_norm(:,didx)); 
%     dmat_lc_norm(:,didx) = dmat_lc_norm(:,didx) / sum(dmat_lc_norm(:,didx)); 
% %     dmat_norm(:,didx) = dmat_norm(:,didx) - min(dmat_norm(:,didx)); 
% %     dmat_norm(:,didx) = dmat_norm(:,didx) / max(dmat_norm(:,didx)); 
% %     dmat_lc_norm(:,didx) = dmat_lc_norm(:,didx) - min(dmat_lc_norm(:,didx)); 
% %     dmat_lc_norm(:,didx) = dmat_lc_norm(:,didx) / max(dmat_lc_norm(:,didx)); 
% 
% end 
% 
% 
% nan(0); 


end

