

tc = toCluster(filter_wasClustered,:); 
tc_mag2 = sum(tc.^2,2); 
[vtc,dtc] = eig(tc'*tc); 
[dtc,sortfilter] = sort(real(diag(dtc)),'descend'); 
vtc = vtc(:,sortfilter); 
ptc = vtc'*tc'; 
filter_far = tc_mag2 > median(tc_mag2); 

% Often get better results by projecting onto a sphere before taking mean 
ptc_sphere = ptc ./ repmat(sqrt(tc_mag2)',[size(ptc,1) 1]); 
ptc_smean = nanmean(ptc_sphere,2); 
ptc_l2 = sqrt(sum((ptc_sphere-repmat(ptc_smean,[1 size(ptc_sphere,2)])).^2,1)); 
ptc_smean = ptc_smean / sqrt(sum(ptc_smean.^2)); 
ptc_cosSphere = ptc_smean'*ptc_sphere; 
ptc_cosSphere(isnan(ptc_cosSphere)) = 0; 

figure(7) 
clf 
hold all 
plot3(ptc(3, filter_far),ptc(2, filter_far),ptc(1, filter_far),'.'); 
plot3(ptc(3,~filter_far),ptc(2,~filter_far),ptc(1,~filter_far),'.'); 
% axis equal 
% axis vis3d 

% Find pcs of "most shaped" (furthest) points from cluster 
filter_far = tc_mag2 > median(tc_mag2); 
ptc_far = ptc(:,filter_far); 
ptc_far_mean = mean(ptc_far,2); 
ptc_far = ptc_far - repmat(ptc_far_mean,[1 size(ptc_far,2)]); 
ptc_centeredAtFar = ptc - repmat(ptc_far_mean,[1 size(ptc,2)]); 
[vtf,dtf] = eig(ptc_far*ptc_far'); 
[dtf,sortfilter] = sort(real(diag(dtf)),'descend'); 
vtf = vtf(:,sortfilter); 
ptf = vtf'*ptc_centeredAtFar; 

figure(7) 
plot3(ptc_far_mean(3),ptc_far_mean(2),ptc_far_mean(1),'ro') 
plot3(ptc_far_mean(3)+1*[0 vtf(3,1)],ptc_far_mean(2)+1*[0 vtf(2,1)],ptc_far_mean(1)+1*[0 vtf(1,1)],'b-','linewidth',2)
plot3(ptc_far_mean(3)+1*[0 vtf(3,2)],ptc_far_mean(2)+1*[0 vtf(2,2)],ptc_far_mean(1)+1*[0 vtf(1,2)],'r-','linewidth',2)
plot3(ptc_far_mean(3)+1*[0 vtf(3,3)],ptc_far_mean(2)+1*[0 vtf(2,3)],ptc_far_mean(1)+1*[0 vtf(1,3)],'g-','linewidth',2)


% ptf(3,tc_mag<5) = 1; 
figure(8) 
clf 
hold all 
plot3(ptf(1, filter_far),ptf(2, filter_far),ptf(3, filter_far),'.') 
plot3(ptf(1,~filter_far),ptf(2,~filter_far),ptf(3,~filter_far),'.') 
% axis equal 




polyPowers = zeros([3 size(ptf,1)]); 
polyPowers(1 , :) = 0; 
polyPowers(2 , 1) = 1; 
polyPowers(3 , 2) = 1; 
polyPowers(4 , 1) = 2; 
polyPowers(5 , 2) = 2; 
polyPowers(6 , 1) = 1; 
polyPowers(6 , 2) = 1; 
% polyPowers(7 , 1) = 3; 
% polyPowers(8 , 1) = 2; 
% polyPowers(8 , 2) = 1; 
% polyPowers(9 , 1) = 1; 
% polyPowers(9 , 2) = 2; 
% polyPowers(10, 2) = 3; 
srdPolyTerms = srdPolyval(ptf, polyPowers); 
C = eye(size(srdPolyTerms,2)); 
C(2,2) = sum(dtf) / dtf(1); 
C(3,3) = sum(dtf) / dtf(2); 
C(4,4) = (sum(dtf)*sum(dtf)) / (dtf(1)*dtf(1));
C(5,5) = (sum(dtf)*sum(dtf)) / (dtf(2)*dtf(2));
C(6,6) = (sum(dtf)*sum(dtf)) / (dtf(1)*dtf(2)); 
% C(7,7) = (sum(dtf)*sum(dtf)*sum(dtf)) / (dtf(1)*dtf(1)*dtf(1)); 
% C(4,4) = 1000; 
% C(5,5) = 1000; 
eta = 1; 
lambda = 1e-2; 
w0 = ones([size(srdPolyTerms,2) 1]); 
srdbeq = ptf(3,:); 
% pointWeights = (1-exp(-(tc_mag2)/max(tc_mag2)));
% pointWeights = max(ptc_cosSphere,0)'.^1; 
% pointWeights(ptf(3,:)<-0.2) = 0; 
% pointWeights(pointWeights<(1-findKnee(1-pointWeights))) = 0; 
pointWeights = exp(-5*ptc_l2/max(ptc_l2))';
pointWeights(isnan(pointWeights)) = 0;
pointWeights(pointWeights<median(pointWeights)) = 0;
% [xk, wk] = srdsr3( srdPolyTerms(filter_far,:), C, srdbeq(filter_far)', pointWeights, eta, lambda, w0, 1e-6, 1e4 ); 
[xk, wk] = srdsr3_mod( srdPolyTerms, C, srdbeq', pointWeights, eta, lambda, w0, 1e-6, 1e4 ); 
bRetrieved = (srdPolyTerms*xk)'; 


spt_norm = srdPolyTerms./repmat(sqrt(sum(srdPolyTerms.^2,1)),[size(srdPolyTerms,1) 1]);
srdbeq_norm = srdbeq/sqrt(sum(srdbeq.^2));
temp = [spt_norm,srdbeq_norm'];

pp_maxDim = 5; 
pp_pureDimIndices = nan([pp_maxDim 1]); 
pppdidx= 1; 
pp = zeros([1+2*pp_maxDim+pp_maxDim*(pp_maxDim-1)/2, size(ptf,1)]); 
ppidx = 1; 
ppidx = ppidx+1; 
for pidx1 = 1:pp_maxDim 
    pp(ppidx,pidx1) = pp(ppidx,pidx1) + 1; 
    pp_pureDimIndices(pppdidx) = ppidx; 
    pppdidx = pppdidx+1; 
    ppidx = ppidx+1; 
    for pidx2 = pidx1:pp_maxDim 
        pp(ppidx,pidx1) = pp(ppidx,pidx1) + 1;
        pp(ppidx,pidx2) = pp(ppidx,pidx2) + 1;
        ppidx = ppidx + 1; 
    end 
end
candidateTerms = srdPolyval(ptf,pp);
for pidx = 1:numel(pointWeights); 
    candidateTerms(pidx,:) = pointWeights(pidx) * candidateTerms(pidx,:); 
end 
candidateTerms = candidateTerms ./ repmat(sqrt(sum(candidateTerms.^2,1)),[size(candidateTerms,1) 1]); 
candidateTerms_cor = abs(candidateTerms'*candidateTerms); 
% Kill terms that divide possible candidate terms 
for pidx1 = 1:size(candidateTerms_cor,1) 
    hasTermsInCommon = any(repmat(pp(pidx1,:)>0,[size(pp,1) 1]) & pp>0, 2); 
    candidateTerms_cor(pidx1,hasTermsInCommon) = 0; 
end 
candidateTerms_cor(1,:) = 0; 
candidateTerms_cor(:,1) = 0; 


figure(9) 
residual = abs(ptf(3,:)-bRetrieved); 
scatter3(ptf(1,:),ptf(2,:),ptf(3,:),2,pointWeights) 
scatter3(ptf(1,:),ptf(2,:),ptf(3,:),2,exp(-5*residual/max(residual))) 

pointWeights = exp(-5*residual/max(residual))'; 
[xk2, wk2] = srdsr3_mod( srdPolyTerms, C, srdbeq', pointWeights, eta, lambda, w0, 1e-6, 1e4 ); 
bRetrieved2 = (srdPolyTerms*xk2)'; 
figure(8) 
plot3(ptf(1,:),ptf(2,:),bRetrieved2,'o') 

xl = xlim; 
yl = ylim; 
[xg,yg] = meshgrid(linspace(xl(1),xl(2),20),linspace(yl(1),yl(2),20)); 
gridvals = reshape(srdPolyval([xg(:)';yg(:)'],polyPowers)*xk2, size(xg)); 
surf(xg,yg,gridvals,'facecolor',[0.5 0.2 0.2],'facealpha',0.3)


