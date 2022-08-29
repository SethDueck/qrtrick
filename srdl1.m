


% % Try to understand matlab's linprog: 
% % y = m*x+b 
% x = linspace(0,1,201)'; 
% y = 5*x + 17; 
% y = y + 1*randn(size(y)); % Add noise 
% y(end) = -10; % Add outlier 
% % Do linear programming 
% srdf   = [ 0; 0; ones(size(x)); ones(size(x)) ]; % x*0, b*0, posResidual*1, negResidual*1
% srdlb  = [ -inf(1); -inf(1); zeros(size(x)); zeros(size(x)) ]; 
% srdAeq = [ x, ones(size(x)), eye(numel(x)), -eye(numel(x)) ]; % x*optsol(1), 1*optsol(2), +posResidual, -negResidual 
% srdbeq = y'; 
% srdx0 = []; 
% srdopt = optimset('Display','final'); 
% tic 
% optsol = linprog(srdf,[],[],srdAeq,srdbeq,srdlb,[],srdx0,srdopt); 
% toc 


rawL1Data = pThisCentroid; 
% rawL1Data = ptc; 

fitWeightFn = 1 - (distFromP(pointsInThisCluster).^2)'/max(distFromP(pointsInThisCluster).^2); 

% temp = (rawL1Data(end-0,:)-mean(rawL1Data(end-0,:))).^2 + (rawL1Data(end-1,:)-mean(rawL1Data(end-1,:))).^2; 
% fitWeightFn = 1 - temp'/max(temp); 

% Choose number of dimensions to rotate into this-centroid PCs 
numDimsThisFit = 3; 
centeredPThisCentroid = rawL1Data(end-numDimsThisFit+1:end,:); 
centeredPThisCentroid = centeredPThisCentroid - repmat(mean(centeredPThisCentroid,2), [1 size(centeredPThisCentroid,2)]); 
[vrl1d,drl1d] = eig(centeredPThisCentroid*centeredPThisCentroid'); 
rotatedPThisCentroid = vrl1d'*centeredPThisCentroid;

riThisCentroid = sort(randperm(size(rotatedPThisCentroid,2),min(1000,size(rotatedPThisCentroid,2)))); % Take N random datapoints 
ridThisCentroid = squareform(pdist(rotatedPThisCentroid(:,riThisCentroid)')); 
ridthisCentroid_numPixelsToKeep = 30; 
ridthisCentroid_bestPixels = srd_findMostInformativePixels(ridThisCentroid,ridthisCentroid_numPixelsToKeep); 
refPixelsThisCentroid = riThisCentroid(ridthisCentroid_bestPixels); 
distToRefThisCentroid = pdist2(rotatedPThisCentroid',rotatedPThisCentroid(:,refPixelsThisCentroid)');
dtrtcn = distToRefThisCentroid; 
% dtrtcn = dtrtcn - repmat(mean(dtrtcn,1),[size(dtrtcn,1) 1]); 
dtrtcn = dtrtcn - repmat(mean(dtrtcn,2),[1 size(dtrtcn,2)]); 
[vdtrtcn,ddtrtcn] = eig(dtrtcn'*dtrtcn); 
pdtrtcn = vdtrtcn'*dtrtcn'; 
% Detect subgroups 
% kThisCentroid = kmeans(pdtrtcn(end,:)',2); 

% Nearest neighbour search, reuse random indices from QR trick step 
% nnSearchIndices = 1:size(rotatedPThisCentroid,2); % Use all datapoints 
nnSearchIndices   = 1:numel(riThisCentroid); 
% nnSearch_i        = randperm(size(rotatedPThisCentroid,2),1); % Use all datapoints
% nnSearch_i        = randperm(numel(nnSearchIndices),1); % Choose random elements of subsample from QR trick 
[~,nnSearch_i] = min(sum(ridThisCentroid)); % Choose, in a naiive sense, the most well-neighboured point 
nnSearchThisIndex = nnSearch_i; 
listNearestDistances = nan([numel(nnSearchIndices) 1]);
listNearestDistances(1) = 0; 
listNearestPoints = nan([size(nnSearchIndices) 1]); 
listNearestPoints(1) = nnSearch_i; 
nnSearchDistances = inf(size(nnSearchIndices)); 
% distsThisCentroid should be equivalent to ridThisCentroid (without a
%     squareform operation) -- I don't know yet whether using the ri subset
%     works, so I avoid the squareform operation and optimize to find the
%     triu indices from pdist in case this part of the code requires the
%     points to be dense. 
% distsThisCentroid = pdist(rotatedPThisCentroid'); 
distsThisCentroid = pdist(rotatedPThisCentroid(:,riThisCentroid)'); 

% figure(2) 
% clf 
% subplot(4,1,1:3) 
% hold all 
% plot3(rotatedPThisCentroid(end-2,:),rotatedPThisCentroid(end-1,:),rotatedPThisCentroid(end-0,:),'.')
% axis equal 
% subplot(4,1,4) 
% plot(listNearestDistances) 
% xlim([0 numel(listNearestDistances)])
nnSearch_onesVec = ones([numel(nnSearchIndices) 1]); 
nnSearch_neg1Vec = [0 cumsum((numel(nnSearchIndices)-2):-1:0)]; 
for pidx = 1:numel(listNearestDistances)-1 

% subplot(4,1,1:3) 
% plot3(rotatedPThisCentroid(end-2,nnSearchThisIndex),rotatedPThisCentroid(end-1,nnSearchThisIndex),rotatedPThisCentroid(end-0,nnSearchThisIndex),'ro')
% title(num2str(pidx))     
% subplot(4,1,4) 
% plot(listNearestDistances) 
% xlim([0 numel(listNearestDistances)])

    nnSearchIndices(nnSearch_i) = []; 
    nnSearchDistances(nnSearch_i) = []; 
    
    % Method 1: Find distances from original pdist array 
    triangularIndices0 = nnSearch_neg1Vec(1:(nnSearchThisIndex-1)) + (nnSearchThisIndex-1); 
    triangularIndices1 = (nnSearch_neg1Vec(nnSearchThisIndex) + nnSearchThisIndex) + (0:(size(rotatedPThisCentroid,2)-nnSearchThisIndex-1));
    triangularIndices_full = [triangularIndices0 0 triangularIndices1]; 
    triangularIndices = triangularIndices_full(nnSearchIndices); 
    nnSearchDistances_new = distsThisCentroid(triangularIndices); 
    
%     % Method 2: Find distances on the fly using pdist2 
%     pToSearch = rotatedPThisCentroid(:,nnSearchIndices)'; 
%     nnSearchDistances_new = pdist2(pToSearch,rotatedPThisCentroid(:,nnSearchThisIndex)');

    % Running 
    nnSearchDistances = min(nnSearchDistances,nnSearchDistances_new); 
    [nnSearch_m,nnSearch_i] = min(nnSearchDistances); 
    nnSearchThisIndex = nnSearchIndices(nnSearch_i); 
    listNearestDistances(pidx+1) = nnSearch_m; 
    listNearestPoints(   pidx+1) = nnSearchThisIndex; 
end 

% Identify spikes in nearest distance
% ldnKneeFnData = sort(listNearestDistances); 
lndKneeFnData = listNearestDistances; 
lnd_runningMean1LR =  cumsum(       lndKneeFnData'   )  ./ (1:numel(lndKneeFnData)); 
lnd_runningMean2LR =  cumsum(       lndKneeFnData'.^2)  ./ (1:numel(lndKneeFnData)); 
lnd_runningMean1RL = fliplr( ...
                            cumsum(fliplr(lndKneeFnData'   )) ./ (1:numel(lndKneeFnData)) ... 
                            ); 
lnd_runningMean2RL = fliplr( ... 
                            cumsum(fliplr(lndKneeFnData'.^2)) ./ (1:numel(lndKneeFnData)) ... 
                            ); 
lndKneeFn =  sqrt( lnd_runningMean2LR-lnd_runningMean1LR.^2 ) ...
       + sqrt( lnd_runningMean2RL-lnd_runningMean1RL.^2 ); 
   

% Create polynomial for fitting surface 
% polyPowers = zeros([6 size(rawL1Data,1)]); 
polyPowers = zeros([6 size(rotatedPThisCentroid,1)]); 
polyPowers(1 , end-0) = 1; 
polyPowers(2 , end-1) = 1; 
polyPowers(3 ,     :) = 0; 
polyPowers(4 , end-0) = 2; 
polyPowers(5 , end-1) = 2; 
polyPowers(6 , end-0) = 1; 
polyPowers(6 , end-1) = 1; 
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

% srdPolyTerms = srdPolyval(rawL1Data, polyPowers); 
srdPolyTerms = srdPolyval(rotatedPThisCentroid, polyPowers); 

% srdtarget = rawL1Data(end-2,:);
srdtarget = rotatedPThisCentroid(end-2,:);

maxNumElementsInL1 = 800; 
if( maxNumElementsInL1 < size(srdPolyTerms,1) ) 
    % Use up to maxNumElementsInL1 -- take random sample but maintain order
    filter_selectRandomUnique = sort(randperm(size(srdPolyTerms,1), maxNumElementsInL1)); 
    x      = srdPolyTerms(filter_selectRandomUnique,:); 
    srdbeq = srdtarget(filter_selectRandomUnique); 
    f_residualWeights = fitWeightFn(filter_selectRandomUnique); 
else 
    % Use whole dataset 
    x      = srdPolyTerms; 
    srdbeq = srdtarget; 
    f_residualWeights = fitWeightFn; 
end

% srdf = [ zeros([size(x,2) 1]); ones([size(x,1) 1]); ones([size(x,1) 1]) ]; % p_1*0, p_2*0, p_1*p_2*0, 1*0, posResidual*1, negResidual*1
srdf = [ zeros([size(x,2) 1]); f_residualWeights; f_residualWeights ]; % p_1*0, p_2*0, p_1*p_2*0, 1*0, posResidual*1, negResidual*1
srdlb = [ -inf([size(x,2) 1]); zeros([size(x,1) 1]); zeros([size(x,1) 1]) ]; 
srdAeq = [ x, eye(size(x,1)), -eye(size(x,1)) ]; % x*optsol, +posResidual, -negResidual 
srdx0 = []; 
srdopt = optimset('Display','off');
% srdopt = optimset('Display','iter'); 

% % Do fit using built-in L1 method 
% optsol = linprog(srdf,[],[],srdAeq,srdbeq,srdlb,[],srdx0,srdopt); 
% bRetrieved = (srdPolyTerms*optsol(1:size(x,2)))'; 
% residuals = abs(srdtarget-bRetrieved); 

C = eye(size(x,2)); 
eta = 1; 
lambda = 1e-4; 
lambda = 1e-11; 
w0 = ones([size(x,2) 1]); 
[xk, wk] = srdsr3( x, C, srdbeq', eta, lambda, w0, 1e-6, 1e4 ); 
% optsol(1:numel(xk)) = xk; 
bRetrieved = (srdPolyTerms*xk)'; 
residuals = abs(srdtarget-bRetrieved); 

% Find knee of fit residuals -- partition into set that maximizes 
%  std(subset1)+std(subset2); 
residuals_sort = sort(residuals); 
residuals_runningMean1LR =  cumsum(       residuals_sort   )  ./ (1:numel(residuals_sort)); 
residuals_runningMean2LR =  cumsum(       residuals_sort.^2)  ./ (1:numel(residuals_sort)); 
residuals_runningMean1RL = fliplr( ...
                            cumsum(fliplr(residuals_sort   )) ./ (1:numel(residuals_sort)) ... 
                            ); 
residuals_runningMean2RL = fliplr( ... 
                            cumsum(fliplr(residuals_sort.^2)) ./ (1:numel(residuals_sort)) ... 
                            ); 
residualKneeFn =  sqrt( residuals_runningMean2LR-residuals_runningMean1LR.^2 ) ...
       + sqrt( residuals_runningMean2RL-residuals_runningMean1RL.^2 ); 
% [~, residualKneeIdx] = max(residualKneeFn(1:end-2)); 
residualKneeIdx = find(0>diff(residualKneeFn),1,'first'); 
residualKneeVal = residuals_sort(residualKneeIdx); 
residualMaxDifference = 1e-1; 
filter_fitsWell = (residualMaxDifference>abs(residuals)) | (residuals <  residualKneeVal); 
filter_fitsPoor = ~filter_fitsWell; 

if( 0 < numel(xk) )
figure(8) 
clf 
subplot(8,1,1:3) 
hold all 
scatter3(rotatedPThisCentroid(end-0,filter_fitsWell),rotatedPThisCentroid(end-1,filter_fitsWell),bRetrieved(filter_fitsWell),15,fitWeightFn(filter_fitsWell),'o','MarkerFaceColor','g','MarkerFaceAlpha',0.4,'LineWidth',2) 
scatter3(rotatedPThisCentroid(end-0,filter_fitsPoor),rotatedPThisCentroid(end-1,filter_fitsPoor),bRetrieved(filter_fitsPoor),15,fitWeightFn(filter_fitsPoor),'o','MarkerFaceColor','r','MarkerFaceAlpha',0.4,'LineWidth',2) 
scatter3(rotatedPThisCentroid(end-0,filter_fitsWell),rotatedPThisCentroid(end-1,filter_fitsWell),srdtarget(filter_fitsWell),35,'k.','MarkerEdgeAlpha',0.6) 
scatter3(rotatedPThisCentroid(end-0,filter_fitsPoor),rotatedPThisCentroid(end-1,filter_fitsPoor),srdtarget(filter_fitsPoor),35,'r.','MarkerEdgeAlpha',0.6) 
title([num2str(kidx) ': ' num2str(size(x,1)) '  (' num2str(numel(bRetrieved)) ', captures ' num2str(sum(filter_fitsWell)) ')']) 
axis equal 
xl = xlim; 
yl = ylim; 
zl = zlim; 
numGridPoints = 20; 
[XL,YL] = meshgrid(linspace(xl(1),xl(2),numGridPoints),linspace(yl(1),yl(2),numGridPoints)); 
xg = srdPolyval([XL(:),YL(:)]', polyPowers(:,end:-1:end-1));
% zplane = XL*optsol(1) + YL*optsol(2) + 1*XL.*YL*optsol(3) + optsol(4); 
zplane = reshape(xg*xk,size(XL)); 
clear hlist
hlist = nan(2,1); 
surf(XL,YL,zplane,ones(size(XL)),'FaceColor',[0.3 0.3 1],'FaceAlpha',0.3, 'EdgeColor','none'); 
xlim(xl) 
ylim(yl) 
zlim(zl) 
hlist(1) = gca; 

subplot(8,1,4:6) 
hold all 

scatter3(rotatedPThisCentroid(end-0,filter_fitsWell),rotatedPThisCentroid(end-1,filter_fitsWell),bRetrieved(filter_fitsWell)-srdtarget(filter_fitsWell),15,fitWeightFn(filter_fitsWell),'o','MarkerFaceColor','g','MarkerFaceAlpha',0.4,'LineWidth',0.5) 
scatter3(rotatedPThisCentroid(end-0,filter_fitsPoor),rotatedPThisCentroid(end-1,filter_fitsPoor),bRetrieved(filter_fitsPoor)-srdtarget(filter_fitsPoor),15,fitWeightFn(filter_fitsPoor),'o','MarkerFaceColor','r','MarkerFaceAlpha',0.4,'LineWidth',0.5) 
surf(XL,YL,zeros(size(XL)),'FaceColor',[0.3 0.3 1],'FaceAlpha',0.3, 'EdgeColor','none'); 
axis equal 
xlim(xl) 
ylim(yl) 
hlist(2) = gca; 

hlink = linkprop(hlist,{'CameraPosition','CameraUpVector'}); 
rotate3d on 
view(45,10); 

subplot(8,1,7) 
hold all 
yyaxis left
plot(listNearestDistances)
yyaxis right
plot((lndKneeFn)) 

subplot(8,1,8) 
hold all 
plot((residuals_sort-min(residuals_sort))/(max(residuals_sort)-min(residuals_sort))); 
plot((residualKneeFn-min(residualKneeFn))/(max(residualKneeFn)-min(residualKneeFn))); 


end 

