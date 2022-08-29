

% candidateTerms = srdPolyval(ptf,pp); 
% for pidx = 1:numel(pointWeights); 
%     candidateTerms(pidx,:) = pointWeights(pidx) * candidateTerms(pidx,:); 
% end 

% lambda = 1e-2; 
% [xk, wk] = srdsr3_mod( srdPolyTerms, C, srdbeq', pointWeights, eta, lambda, w0, 1e-6, 1e4 );

% temp0 = ptc_sphere(1:3,:)'; 
% for tidx = 1:1000
%     temp1 = temp0; 
%     for didx = 1:3 
%         temp1(:,didx) = temp1(:,didx)-ptc_sphere(didx,10); 
%     end 
%     temp2 = sum(temp1.^2,2); 
% %     temp2 = sum((ptc_sphere(1:3,:) - repmat(ptc_sphere(1:3,10),[1 size(ptc_sphere,2)])).^2); 
% end 

% [u,s,v] = srd_svd(ssdata,0e-3);

if(~exist('mdl','var')) 
    mdl = KDTreeSearcher(ssdata,'BucketSize',10);
end
if(~exist('nnIdx','var')) 
%     [nnIdx,nnD] = knnsearch(mdl,ssdata(:,:),'k',20);
    [nnIdx,nnD] = knnsearch(mdl,ssdata(:,:),'k',50);
end 
% b = KnnFind.Approximate(ssdata,ssdata,'K',20,'NSMethod','nn_descent');

[vs,ds] = eig(ssdata'*ssdata); 
[ds,sortfilter] = sort(real(diag(ds)),'descend'); 
vs = vs(:,sortfilter); 
vs = vs(:,1:20); 
ps = ssdata*vs; 
if(~exist('mdl2','var')) 
    mdl2 = KDTreeSearcher(ps,'BucketSize',10);
end
if(~exist('nnIdx2','var')) 
    [nnIdx2,nnD2] = knnsearch(mdl2,ps,'k',20);
end 

rp_axes = 2*(0.5-rand([10,size(ssdata,2)])); 
for ridx = 1:size(rp_axes,1) 
    rp_axes(ridx,:) = rp_axes(ridx,:) / sqrt(sum(rp_axes(ridx,:).^2)); 
end 
rp = ssdata*rp_axes'; 

needsCluster = true([size(ssdata,1) 1]); 
clusterList = zeros([size(ssdata,1) 1]); 
clusterCovariance = zeros(3*ceil(size(ssdata,1)/size(nnIdx,2)) * [1 1]); 
hitCount1 = 0; 
hitCount2 = 0; 
clusterIdx = 0; 
[~,goodNeighbourOrder] = sort(sum(nnD,2)); % Start with the points with the closest neighbours 
for cidx_linear = 1:numel(needsCluster) 
    cidx = goodNeighbourOrder(cidx_linear); 
    if(needsCluster(cidx)) 
        clusterIdx = clusterIdx + 1; 
        for nidx = 1:size(nnIdx,2) 
            if(needsCluster(nnIdx(cidx,nidx))) 
                needsCluster(nnIdx(cidx,nidx)) = false; 
                clusterList(nnIdx(cidx,nidx)) = clusterIdx; 
                hitCount1 = hitCount1+1; 
            else 
                clusterCovariance(clusterList(nnIdx(cidx,nidx)),clusterIdx) = ...
                clusterCovariance(clusterList(nnIdx(cidx,nidx)),clusterIdx) + 1; 
                hitCount2 = hitCount2 + 1; 
            end
        end 
    end 
end 
clusterCovariance = clusterCovariance + clusterCovariance'; 
filter = 0==sum(clusterCovariance);
clusterCovariance(filter,:) = []; 
clusterCovariance(:,filter) = []; 
temp = clusterCovariance; 
for tidx = 1:size(temp,2)
    temp(:,tidx) = temp(:,tidx)/sum(temp(:,tidx));
end
[v,d] = eig(temp'*temp); 
[d,sortfilter] = sort(diag(d),'descend'); 
v = v(:,sortfilter); 
p = v'*temp';
numDimsToSearch = find(d>1,1,'last'); 
[m1,clist] = max(abs(p(1:numDimsToSearch,:))); 
[clist_unique,~,clist_reduced] = unique(clist);

clabels = nan([size(ssdata,1) 1]); 
for cidx = 1:numel(clist_unique)
	clabels(ismember(clusterList,find(cidx==clist_reduced))) = cidx; 
end 

% Make sure clusters have good shape 
methodsToTry.gaps   = true; 
methodsToTry.tails  = true; 
methodsToTry.unfold = true; 
cidx = 31; 
while cidx <= max(clabels)
    % Test if cluster has good shape 
    thisCluster = ssdata(clabels==cidx,:); 
    thisCluster = thisCluster - mean(thisCluster(:)); 
    [vc,dc] = eig(thisCluster'*thisCluster); 
    [dc,sortfilter] = sort(real(diag(dc)),'descend'); 
    vc = vc(:,sortfilter); 
    pc = vc(:,(dc/sum(dc))>1e-3)'*thisCluster'; 
%     pc = vc(:,dc>1)'*thisCluster'; 
    debugFlag = false; 
    if( debugFlag || ( 1-((dc(1)+dc(2))/sum(dc))) > 1e-2 ) 
        % Cluster has bad shape, needs to be split up 
        [subclusterIndices,methodsToTry] = nnCluster_splitCluster( pc, thisCluster, methodsToTry ); 
        if( 1 < max(subclusterIndices) )
            % Cluster can be split. We'll relabel elements of the cluster
            % and try again. 
            clabels = renameClusterIndices(clabels, subclusterIndices, cidx); 
            continue; 
        else 
            % Couldn't find a way to split cluster -- no choice but to
            % continue without improving cluster uniformity 
            warning('Had to continue without splitting cluster'); 
        end 
    end 
        
    
    
    % Test if cluster has outliers 
    thisCluster = thisCluster - repmat(mean(thisCluster,1),[size(thisCluster,1) 1]);
    [vcn,dcn] = eig(thisCluster'*thisCluster); 
    [dcn,sortfilter] = sort(real(diag(dcn)),'descend'); 
    vcn = vcn(:,sortfilter); 
    
    % Move on to next cluster 
    cidx = cidx + 1; 
    methodsToTry.gaps   = true; 
    methodsToTry.tails  = true; 
    methodsToTry.unfold = true; 
    
end 

% From wed night: 
ptf_n = ptf(1:5,:); 
for didx = 1:size(ptf_n,1) 
    ptf_n(didx,:) = ptf_n(didx,:) / std(ptf_n(didx,:)); 
end 
[k,c] = kmeans(ptf_n(1:5,:)',20);
t = delaunay([c(:,1);0],[c(:,2);0],[c(:,3);0]);
su = turnTriIntoSurface( t ); 
sunz = su; 
sunz(any(sunz>size(c,1),2),:) = [];

sue = turnTriIntoSurface_edges( t ); 
suenz = sue; 
suenz(any(suenz>size(c,1),2),:) = [];

edgeScores = nan([1 1]*max(k)); 
for kidx1 = 1:max(k) 
    edgeScores(kidx1,kidx1) = 0; 
    for kidx2 = kidx1+1:max(k) 
dnn = findNearestNeighbour(ptc_sphere(1:3,k==kidx1)',ptc_sphere(1:3,k==kidx2)'); 
edgeScores(kidx1,kidx2) = dnn; 
edgeScores(kidx2,kidx1) = dnn; 
    end 
end 
edgeScores_reduced = reduceConnectivityGraph(edgeScores); 

es_n = edgeScores; 
esr_n = edgeScores_reduced; 
for eidx = 1:size(esr_n,2) 
    es_n(:,eidx) = exp(-5*es_n(:,eidx)/max(es_n(:,eidx))); 
    esr_n(:,eidx) = exp(-5*esr_n(:,eidx)/max(esr_n(:,eidx))); 
%     esr_n(:,eidx) = esr_n(:,eidx) / sum(esr_n(:,eidx)); 
end 
es_n = min(es_n,es_n'); 
esr_n = min(esr_n,esr_n'); 

figure(411) 
clf 
hold all 
for kidx = 1:max(k) 
    plot3(ptf_n(1,kidx==k),ptf_n(2,kidx==k),ptf_n(3,kidx==k),'.')
end 
faceScores = nan([size(sunz,1) 1]); 
edgeScoresDelta = edgeScores - edgeScores_reduced; 
% edgeScoresDelta = 1 - (edgeScoresDelta/max(edgeScoresDelta(:))); 
for tidx = 1:size(sunz,1); 
%     faceScores(tidx) = edgeScoresDelta(sunz(tidx,1),sunz(tidx,2)) * ... 
%                        edgeScoresDelta(sunz(tidx,1),sunz(tidx,3)) * ... 
%                        edgeScoresDelta(sunz(tidx,2),sunz(tidx,3)); 
    
    faceScores(tidx) = esr_n(sunz(tidx,1),sunz(tidx,2)) * ... 
                       esr_n(sunz(tidx,1),sunz(tidx,3)) * ... 
                       esr_n(sunz(tidx,2),sunz(tidx,3)); 
%    faceScores(tidx) = min([esr_n(sunz(tidx,1),sunz(tidx,2)), ... 
%                        esr_n(sunz(tidx,1),sunz(tidx,3)), ... 
%                        esr_n(sunz(tidx,2),sunz(tidx,3))]); 
end 
faceScores =   faceScores/max(faceScores); 
for tidx = 1:size(sunz,1)
    trimesh(sunz(tidx,:),[c(:,1);0],[c(:,2);0],[c(:,3);0],'edgecolor','k','facecolor',[1 1 1]*0.2,'facealpha',faceScores(tidx)*0.4+0 )
end 
% sueScore = nan([size(suenz,1) 1]); 
% for eidx = 1:size(suenz,1) 
%     sueScore(eidx) = edgeScoresDelta(suenz(eidx,1),suenz(eidx,2)); 
% end 
% sueScore =  sueScore/max(sueScore); 
% for eidx = 1:size(suenz,1) 
%     plot3([c(suenz(eidx,1),1),c(suenz(eidx,2),1)],[c(suenz(eidx,1),2),c(suenz(eidx,2),2)],[c(suenz(eidx,1),3),c(suenz(eidx,2),3)],'color',[1 1 1]*sueScore(eidx)*0.2+0.8); 
% end 


% 
% [k,c] = kmeans(ptc_sphere(1:5,:)',20);
% edgeScores = nan([1 1]*max(k)); 
% for kidx1 = 1:max(k) 
%     edgeScores(kidx1,kidx1) = 0; 
%     for kidx2 = kidx1+1:max(k) 
% dnn = findNearestNeighbour(ptc_sphere(1:3,k==kidx1)',ptc_sphere(1:3,k==kidx2)'); 
% edgeScores(kidx1,kidx2) = dnn; 
% edgeScores(kidx2,kidx1) = dnn; 
%     end 
% end 
% edgeScores_reduced = reduceConnectivityGraph(edgeScores); 
c_sphere = nan(size(c)); 
for kidx = 1:max(k) 
    c_sphere(kidx,:) = mean(ptf(1:size(c_sphere,2),kidx==k),2)'; 
end

temp = exp(-5*edgeScores_reduced/max(edgeScores_reduced(:)));
temp2 = temp;
for tidx = 1:size(temp,2)
%     temp2(:,tidx) = temp2(:,tidx)/sum(temp2(:,tidx));
    temp2(:,tidx) = temp2(:,tidx)/sqrt(sum(temp2(:,tidx).^2));
end 
[ves,des] = eig(temp2'*temp2);
pes = ves'*temp2; 
[~,tempidx] = max(abs(pes));
filter_clusterIsInMainBody = mode(tempidx)==tempidx; 
filter_faceIsInMainBody = all(filter_clusterIsInMainBody(sunz),2); 


faceScores_c = nan([size(sunz,1) 1]); 
edgeScoresDelta = edgeScores_reduced ./ edgeScores; 
esDelta = es_n./esr_n; 
% edgeScoresDelta = 1 - (edgeScoresDelta/max(edgeScoresDelta(:))); 
for tidx = 1:size(sunz,1); 
%     faceScores(tidx) = edgeScoresDelta(sunz(tidx,1),sunz(tidx,2)) * ... 
%                        edgeScoresDelta(sunz(tidx,1),sunz(tidx,3)) * ... 
%                        edgeScoresDelta(sunz(tidx,2),sunz(tidx,3)); 
    if( filter_faceIsInMainBody(tidx) ) 
%     faceScores_c(tidx) = esr_n(sunz(tidx,1),sunz(tidx,2)) * ... 
%                        esr_n(sunz(tidx,1),sunz(tidx,3)) * ... 
%                        esr_n(sunz(tidx,2),sunz(tidx,3)); 
    faceScores_c(tidx) = esDelta(sunz(tidx,1),sunz(tidx,2)) * ... 
                       esDelta(sunz(tidx,1),sunz(tidx,3)) * ... 
                       esDelta(sunz(tidx,2),sunz(tidx,3)); 
%     faceScores_c(tidx) = min([esr_n(sunz(tidx,1),sunz(tidx,2)) , ... 
%                        esr_n(sunz(tidx,1),sunz(tidx,3)) , ... 
%                        esr_n(sunz(tidx,2),sunz(tidx,3))]); 
    else 
        faceScores_c(tidx) = 0; 
    end 
    
%    faceScores(tidx) = min([esr_n(sunz(tidx,1),sunz(tidx,2)), ... 
%                        esr_n(sunz(tidx,1),sunz(tidx,3)), ... 
%                        esr_n(sunz(tidx,2),sunz(tidx,3))]); 
end 
faceScores_c(faceScores_c>1) = 1; 

figure(20) 
clf 
hold all 
for kidx = 1:max(k)
    if( filter_clusterIsInMainBody(kidx) )
        plot3(ptf(1,k==kidx),ptf(2,k==kidx),ptf(3,k==kidx),'.')
    else 
%         plot3(ptf(1,k==kidx),ptf(2,k==kidx),ptf(3,k==kidx),'kx') 
        set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')+1); 
    end
end
for tidx = 1:size(sunz,1)
    if( filter_faceIsInMainBody(tidx) )
%         trimesh(sunz(tidx,:),[c_sphere(:,1);0],[c_sphere(:,2);0],[c_sphere(:,3);0],'edgecolor','k','facecolor',[1 1 1]*0.2,'facealpha',0.6*(faceScores(tidx)/max(faceScores(filter_faceIsInMainBody))) )
        trimesh(sunz(tidx,:),[c_sphere(:,1);0],[c_sphere(:,2);0],[c_sphere(:,3);0],'edgecolor','k','facecolor',[1 1 1]*0.2,'facealpha',faceScores_c(tidx) )
%         display(num2str(faceScores(tidx))) 
nan(0); 
    end 
end 
axis normal 


