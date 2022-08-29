

filter = clabels < 108; 

figure(108) 
clf 
hold all 
plot(wavs,log10(hr(:,40))) 
for cidx = 1:98%max(clabels) 
    filter = clabels==cidx; 
    sdcent = ssdata(filter,:) - repmat(mean(ssdata(filter,:)),[sum(filter) 1]); 
    [~,clusterMedIdx] = min(sum(sdcent.^2,2)); 
    sdMedIdx = find(filter,clusterMedIdx,'first'); 
    plot(wavs(filter),log10(hr(sdMedIdx(end),40))*ones([sum(filter) 1]),'.') 
end 



