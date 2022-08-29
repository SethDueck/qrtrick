

% Find subsets that are stable 
dToNeighbour = sqrt(sum( (pdata(:,2:end)-pdata(:,1:end-1)).^2 )); 
dToNeighbour_limit = 0.1; 
filter_fromContinuity = ([0 dToNeighbour] < dToNeighbour_limit) & ([dToNeighbour 0] < dToNeighbour_limit); 
filter_list = find(filter_fromContinuity); 

% Create filter for any datapoints that are also near stable points 
filter_stableNeighbourTestSubset = filter_list(1);  
distToStableSubset_running = pdist2(pdata(:,filter_fromContinuity)',pdata(:,filter_stableNeighbourTestSubset(end))'); 
[mval,midx] = max(distToStableSubset_running); 
while( mval > 10*dToNeighbour_limit ) 
    filter_stableNeighbourTestSubset = [filter_stableNeighbourTestSubset, filter_list(midx)];  %#ok<AGROW>
    distToStableSubset_candidate = pdist2(pdata(end-20:end,filter_fromContinuity)',pdata(end-20:end,filter_stableNeighbourTestSubset(end))'); 
    distToStableSubset_running = min(distToStableSubset_running,distToStableSubset_candidate); 
    [mval,midx] = max(distToStableSubset_running); 
end 
distToStableSubset = pdist2(pdata',pdata(:,filter_stableNeighbourTestSubset)'); 
filter = min(distToStableSubset,[],2) < 10*dToNeighbour_limit; 

data1 = data; 
data1_mean = mean(data(filter,:),1); 
data1 = data1 - repmat(data1_mean,[size(data1,1) 1]); 
[vf,df] = eig(data1(filter,:)'*data1(filter,:)); 
[~,sortfilter] = sort(diag(df)); 
vf = vf(:,sortfilter); 
df = df(sortfilter,sortfilter); 
data0 = data; 
data0 = data0 - repmat(mean(data(~filter,:),1),[size(data0,1) 1]); 
[vnf,dnf] = eig(data0(~filter,:)'*data0(~filter,:)); 
[~,sortfilter] = sort(diag(dnf)); 
vnf = vnf(:,sortfilter); 
dnf = dnf(sortfilter,sortfilter); 

pf11 = vf' *data1( filter,:)'; 
pf10 = vf' *data1(~filter,:)'; 
pf01 = vnf'*data0( filter,:)'; 
pf00 = vnf'*data0(~filter,:)'; 
pf0 = vnf'*data0'; 

% Use QR trick to find best interior points 
interiorDists_orig = squareform(pdist(data(~filter,:)));
interiorDists = interiorDists_orig; 
interiorIdx = nan(10,1); 
for iidx = 1:numel(interiorIdx) 
    interior_mags = sum(interiorDists.^2,2); 
    [~,bestIdx] = max(interior_mags); 
    interiorIdx(iidx) = bestIdx; 
    bestIdx_unit = interiorDists(bestIdx,:); 
    bestIdx_unit = bestIdx_unit / sqrt(sum(bestIdx_unit.^2)); 
    interiorDists = interiorDists- (interiorDists*bestIdx_unit')*bestIdx_unit; 
end 
temp = find(~filter); 
filter_dataToQRSubset = temp(interiorIdx); 

% Take random subset of "interior" points
filter_innerEls = sort(randperm(size(pf00,2),1000)); 
temp = find(~filter); 
filter_dataToInnerEls = temp(filter_innerEls); 
pd = squareform(pdist(pf00(:,filter_innerEls)')); 


[v,d] = eig(pd); 
[~,sortfilter] = sort(abs(diag(d))); 
d = d(sortfilter,sortfilter); 
v = v(:,sortfilter); 
d = d(2:end,2:end); 
v = v(:,2:end); 

temp = pdist2(pf0',pf00(:,filter_innerEls)');
ptemp = v'*temp';

figure(17) 
clf 
hold all 
scatter3(ptemp(end-2,filter),ptemp(end-1,filter),ptemp(end-0,filter)+200,10,hr(filter,30),'s')
hold all
scatter3(ptemp(end-2,~filter),ptemp(end-1,~filter),ptemp(end-0,~filter),5,hr(~filter,30))
plot3(ptemp(end-2,filter_dataToQRSubset),ptemp(end-1,filter_dataToQRSubset),ptemp(end-0,filter_dataToQRSubset),'ko')

% Identify centroids of subsets, compute vector between mean of subsets 
k = kmeans(sum(pf11.^2,1)',2); 
m1 = mean(pf11(:,k==1),2); 
m2 = mean(pf11(:,k==2),2); 

figure(2) 
clf 
subplot(2,2,1) 
scatter3(pdata(end-2, filter),pdata(end-1, filter),pdata(end-0, filter),10,log10(hr( filter,40)))
colorbar
subplot(2,2,2) 
plot(wavs(filter),hr( filter,20),'.') 

subplot(2,2,3) 
scatter3(pdata(end-2,~filter),pdata(end-1,~filter),pdata(end-0,~filter),10,log10(hr(~filter,40)))
colorbar 
subplot(2,2,4) 
plot(wavs(~filter),hr(~filter,20),'.') 

