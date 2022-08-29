


if( ~exist('ke','var') ) 
    load('a-band-data/optprops.mat'); 
end 

data = ks ./ ke; 
% data = profileToCumsumTransform(data,100); 
data = data -  repmat(mean(data,  1),  [size(data,1),1]); 
data = data ./ repmat(std( data,0,1),[size(data,1),1]); 
[vdata,ddata] = eig(data'*data); 

% Choose centroids using lloyd's algorithm 
pdata_nDim = 6; 
pdata = vdata'*data'; 
% pdata = pdata((end-pdata_nDim+1):end,:); 
[k,c] = kmeans(pdata',15,'MaxIter',1000); 

% Choose medoids using greedy QR trick 
randomIndices = sort(randperm(size(pdata,2),1000)); % Take N random datapoints 
riDistances = squareform(pdist(pdata(:,randomIndices)')); 
rid_numPixelsToKeep = 30; 
rid_bestPixels = srd_findMostInformativePixels(riDistances,rid_numPixelsToKeep); 
referenceMedoidIndices = randomIndices(rid_bestPixels); 
dToReference = pdist2(pdata',pdata(:,referenceMedoidIndices)'); 
[~,dtrnk] = min(dToReference,[],2); 
dtrn = dToReference; 
dtrn = dtrn -  repmat(mean(dtrn,  1), [size(dtrn,1),1]); 
dtrn = dtrn ./ repmat(std( dtrn,0,1), [size(dtrn,1),1]); 
[vdtrn,ddtrn] = eig(dtrn'*dtrn); 
pdtrn = vdtrn'*dtrn'; 


% Find distance from each datapoint to each centroid 
kdist = nan([max(k),size(pdata,2)]); 
for kidx = 1:max(k) 
    kdist(kidx,:) = sqrt(sum((pdata-repmat(c(kidx,:)',[1 size(pdata,2)])).^2)); 
end 
kdist = kdist -  repmat(mean(kdist,2),[1 size(kdist,2)]); 
[vdist,ddist] = eig(kdist*kdist');
pkdist = vdist'*kdist;

% kdist = kdist ./ repmat(std(kdist,1,2),[1 size(kdist,2)]); 

% dLibrary = cell([max(k) 1]); 
% for kidx = 1:max(k) 
% % for kidx = 1:200
%     dLibrary{kidx} = squareform(pdist(pdata(:,k==kidx)')); 
% %     dLibrary{kidx} = squareform(pdist(p(:,1:floor(size(p,2)/200))')); 
% end 

% Find mds dirs within each cluster 
clustersToAssign = 1:size(pdata,2); 
% Start with most highly-populated clusters, so I can look at the most
% interesting ones first during debugging 
[~,dtrnk_ordered] = sort(hist(dtrnk,max(dtrnk))); 
dtrnk_ordered = fliplr(dtrnk_ordered); 
for kidx = dtrnk_ordered; 
    
    % Choose points that are in this centers Voronoi region 
%     pThisCentroid = pdata(:,k==kidx); 
    pointsInThisCluster = dtrnk==kidx; 
    pThisCentroid = pdata(:,pointsInThisCluster);
    
% %     Greedy algorithm 
    numPointsInCluster = 1000; 
    randPointIdx = clustersToAssign(randperm(numel(clustersToAssign),1)); 
    distFromP = sum((pdata-repmat(pdata(:,randPointIdx),[1 numel(clustersToAssign)])).^2); 
%     distFromP_kdistSpace = sum((pkdist-repmat(pkdist(:,randPointIdx),[1 numel(clustersToAssign)])).^2); 
%     distFromP_qrSpace = sum( (pdtrn-repmat(pdtrn(:,randPointIdx),[1 numel(clustersToAssign)])).^2); 
%     
% %     [~,pointDistanceOrder] = sort(distFromP); 
%     [~,pointDistanceOrder] = sort(distFromP_qrSpace); 
%     pointsInThisCluster = pointDistanceOrder(1:numPointsInCluster); 
%     pThisCentroid = pdata(:,pointsInThisCluster); 
%     distThisCentroid = squareform(pdist(pThisCentroid')); 
    
%     % Find distances between points in this set 
%     m = squareform(pdist(pThisCentroid')); 
%     
%     % Double centering 
% %     A = -m/2; 
% %     n = size(pThisCentroid,2); 
% %     C = eye(n) - (1/n)*ones([1 1]*n); 
% %     B_orig = (-1/2)*C*m*C; 
%     B = (-1/2)*m; 
%     B_mean = mean(B,2); 
%     for bidx = 1:size(B,2) 
%         B(:,bidx) = B(:,bidx) - B_mean(bidx);
%     end 
%     B = B';
%     B_mean = B_mean - mean(B_mean); 
%     for bidx = 1:size(B,2) 
%         B(:,bidx) = B(:,bidx) - B_mean(bidx);
%     end 
%     
%     [v,d] = eigs(B,5); 
%     [~,filterdimorder] = sort(diag(real(d))); 
%     v = v(:,filterdimorder); 
%     d = d(filterdimorder,filterdimorder); 
%     ptc = v'*B; 
    srdl1; 
    dimOutlierRank = sum(abs(v).^4); 
    pOutlierRank = dimOutlierRank*(ptc.^1)';
    nan(0); 
end 

m = nan([size(p,2),size(c,1)]); 
for didx = 1:size(p,2) 
    temp = c - repmat(p(:,didx)',[size(c,1),1]); 
    m(didx,:) = sqrt(sum(temp.^2,2)); 
end
[vm,dm] = eig(m'*m); 
pm = vm'*m'; 

m_norm = m; 
m_norm = m_norm - repmat(mean(m_norm,1),[size(m_norm,1) 1]); 
m_norm = m_norm ./ repmat(std(m_norm,0,1),[size(m_norm,1) 1]); 
% m_norm = m_norm - repmat(mean(m_norm,2),[1 size(m_norm,2)]); 
% m_norm = m_norm ./ repmat(std(m_norm,0,2),[1 size(m_norm,2)]); 

m2 = nan([size(c,1),size(c,1)]); 
for didx = 1:size(c,1) 
    temp = c - repmat(c(didx,:),[size(c,1),1]); 
    m2(didx,:) = sqrt(sum(temp.^2,2)); 
end
[vm2,dm2] = eig(m2'*m2); 
pm2 = vm2'*m2'; 
pm_smallOnLarge = vm2'*m';

directions = nan([size(c,1),size(c,1),size(c,2)]); 
for cidx1 = 1:size(c,1) 
    for cidx2 = 1:size(c,1) 
        if( 0 < m2(cidx1,cidx2) )
            temp = c(cidx1,:)-c(cidx2,:); 
            directions(cidx1,cidx2,:) = temp / m2(cidx1,cidx2); 
        else 
            directions(cidx1,cidx2,:) = 0.0; 
        end
    end 
end 


optscalefactor = nan([size(c,1),1]);
numScaleFactors = 100; 
dlist = nan([numScaleFactors,size(c,1)]);
tlist = nan([numScaleFactors,size(c,1)]);
tempr_orig = rand([size(m2,1) 1]); 
tempr_orig = tempr_orig / sqrt(sum(tempr_orig.^2)); 
for cidx = 1:size(c,1) 
   

    
    tempr = tempr_orig; 
    temprm = nan(1); 
    temprm_old = nan(1); 
    distances_thisCenter = m2(:,cidx); 
    smalldirs = reshape(directions(cidx,:,:),[size(directions,2),size(directions,3)]); 
    cosUnits = smalldirs*smalldirs'; 
    
    sfbound_min = min(distances_thisCenter(distances_thisCenter>0));
    sfbound_max = max(distances_thisCenter(distances_thisCenter>0));
    scalefactorlist = logspace(log10(0.001/sfbound_max),log10(5/sfbound_min),100);     
    for sidx = 1:numel(scalefactorlist) 
    
        isNewScalefactor = true; 
%         magsScaled = 10.^(-scalefactorlist(sidx)*distances_thisCenter);
        magsScaled = exp(-scalefactorlist(sidx)*distances_thisCenter);
        magsScaled2 = magsScaled*magsScaled';
        temp2 = magsScaled2.*cosUnits; 
        
%             if(2==cidx || 57==cidx || 73==cidx) 
%                figure(5) 
%                clf 
%                hold all 
%                pt = repmat(magsScaled,[1 3]).*smalldirs(:,end-2:end);
%                plot3(pt(:,1),pt(:,2),pt(:,3),'.');
%                plot3(0,0,0,'ko')
%                axis equal 
%                xlim([-1 1]*1.1)
%                ylim([-1 1]*1.1)
%                zlim([-1 1]*1.1)
%                title(num2str(sidx)) 
%                view(-85,31)
%                drawnow 
%             end
%         [vtemp,dtemp] = eig(temp*temp'); 
%         [vtemp,dtemp] = eigs(temp2,2); % Alternative, maybe faster
%         dlist(sidx) = sum(diag(dtemp)); 
%         dlist(sidx,cidx) = dtemp(end,end)/sum(diag(dtemp)); 
        
        while (isNewScalefactor || ~(abs(abs((temprm/temprm_old))-1) < 1e-6))
            isNewScalefactor = false; 
            temprm_old = temprm; 
            tempr = temp2*tempr; 
            temprm = sqrt(sum(tempr.^2)); 
            if( isfinite(temprm) )
                tempr = tempr / temprm; 
            else 
                display('temprm is nan'); 
                continue
            end
        end 
        dlist(sidx,cidx) = temprm;
        tlist(sidx,cidx) = trace(temp2); 
%         dlist(sidx,cidx) = temprm; 
%         display(num2str(temprm/temprm_old)); 
        
    end
%     [~,midx] = max((diff(dlist(:,cidx)./tlist(:,cidx)))); 
    [~,midx] = max(+diff(dlist(:,cidx)).^2 + diff(tlist(:,cidx).^2 )); 
%     [~,midx] = min(diff((diff(tlist(:,cidx))./diff(scalefactorlist')).^2 - (diff(dlist(:,cidx))./diff(scalefactorlist')).^2));
    optscalefactor(cidx) = scalefactorlist(midx); 
    nan(0); 
end 
    
dtratio = dlist./tlist; 
dtrn = dtratio - repmat(dtratio(1,:),[size(dtratio,1) 1]);
dtrn = dtrn./repmat(sqrt(sum(dtrn.^2,1)),[size(dtrn,1) 1]); 
simmap = dtrn'*dtrn; 
[vdtn,ddtn] = eig(simmap); 
pdtn = vdtn'*dtrn'; 
[~,sortfilter] = sort(pdtn(end,:)); 

figure(2) 
% clf
hold all 
% plot3(pm_smallOnLarge(end-2,:),pm_smallOnLarge(end-1,:),pm_smallOnLarge(end-0,:),'.')
% plot3(pm2(end-2,:),pm2(end-1,:),pm2(end-0,:),'o')
scatter3(pm2(end-2,:),pm2(end-1,:),pm2(end-0,:),20,log(optscalefactor))

