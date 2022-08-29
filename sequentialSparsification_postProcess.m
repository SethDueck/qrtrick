

filterb2 = false([size(ssdata,1) 1]);
filterb2(filter1) = true; 
indexListRound2 = find(~filterb2); 
ssdata2 = ssdata(indexListRound2,:); 

%
% Mostly copied from "first run" code 
%
% ssdata2 = ssdata2 -  repmat(mean(ssdata2,  1), [size(ssdata2,1),1]); 
% ssdata2 = ssdata2 ./ repmat(std( ssdata2,0,1), [size(ssdata2,1),1]); 

[mu2,lambda2] = getDefaultRpcaParams(ssdata2); 

figure(5) 
clf 
subplot(4,1,1)
hold all 
subplot(4,1,2)
hold all 
subplot(4,1,3)
hold all 

errorSequence = 10.^linspace(-2,-5,20); 
% errorSequence = 10.^linspace(-2,-2,30); 
% filter = true([size(ssdata,1) 1]); 
filter2 = 1:size(ssdata2,1); 
filter2_smallSparsePass = true([size(ssdata2,1) 1]); 
S2_old = nan(1);  
Y2_old = nan(1); 
isConverged = false; 
LRank2_old = []; 
for erroridx = 1:numel(errorSequence)
    if( isConverged ) 
        break; 
    end 
    
    filter2 = filter2(filter2_smallSparsePass); 
%     [L1,S1] = rpca(ssdata(filter,:),1e-2, mu, lambda);
    [mu2,lambda2] = getDefaultRpcaParams(ssdata2(filter2,:)); 
    [L2,S2, Y2, LRank2] = rpca(ssdata2(filter2,:),errorSequence(erroridx), mu2, lambda2, S2_old, Y2_old, LRank2_old ); 
    % X = L1 + S1 + Y/mu -- this *may* be true after the first pass 

    unsmoothPart2 = ssdata2(filter2,:) - L2; 
    unsmoothPart2_mag = sqrt(sum(unsmoothPart2.^2,2)); 
    filter2_smallSparsePass = (unsmoothPart2_mag-min(unsmoothPart2_mag)) < ... 
        ((max(unsmoothPart2_mag)-min(unsmoothPart2_mag))/4); 
    
    fracErrorL  = sqrt(sum((ssdata2(filter2,:)-L2   ).^2,2)) ...
        ./ sqrt(sum((ssdata2(filter2,:)).^2,2)); 
    fracErrorLS = sqrt(sum((ssdata2(filter2,:)-L2-S2).^2,2)) ...
        ./ sqrt(sum((ssdata2(filter2,:)).^2,2));     
    isConverged = max(fracErrorLS) < 1e-3; 
    
    
	Y2_old = Y2(filter2_smallSparsePass,:); 
    S2_old = S2(filter2_smallSparsePass,:); 
    L2_old = L2; 
    LRank2_old = LRank2; 
    
    figure(5) 
    subplot(4,1,1)
    plot(filter2, sqrt(sum(S2.^2,2)) ,'.')
    title([num2str(erroridx) ': ' num2str(numel(filter2))]); 
    ylabel('S') 
    subplot(4,1,2) 
    plot(filter2, sqrt(sum((ssdata2(filter2,:)-L2).^2,2)) ,'.')
    ylabel('data - L') 
    nan(0); 
    subplot(4,1,3) 
    plot(filter2, sqrt(sum((ssdata2(filter2,:)-L2-S2).^2,2)) ,'.')
    ylabel('data - L - S') 
    subplot(4,1,4)
    cla 
    hold all 
    plot(filter2, fracErrorL ,'.')
    plot(filter2, fracErrorLS,'.') 
    ylabel('frac error') 
    title(['maxL = ' num2str(max(fracErrorL )) ',   maxLS = ' num2str(max(fracErrorLS))]); 
    
%     figure(5) 
%     clf 
%     plot3(sqrt(sum((ssdata(filter,:)-L1-S1).^2,2)), sqrt(sum((ssdata(filter,:)-L1).^2,2)), sqrt(sum(S1.^2,2)),'.') 
%     xlabel('data - L - S') 
%     ylabel('data - L') 
%     zlabel('S') 
%     axis vis3d

    nan(0); 
end


[u_L2,s_L2,v_L2] = svd(L2,'econ'); 
[u_S2,s_S2,v_S2] = svd(S2,'econ'); 
[u_r2,s_r2,v_r2] = svd(ssdata2(filter2,:)-L2,'econ'); 

p_L2 = v_L2'*L2'; 
p_S2 = v_S2'*S2'; 
p_r2 = v_r2'*(ssdata2(filter2,:)-L2)'; 

hr_subset = hr(~filterb2,:); 
figure(3) 
clf 
subplot(1,2,1) 
scatter3(p_L2(3,:),p_L2(2,:),p_L2(1,:),1,hr_subset(filter2,10))
axis vis3d 
subplot(1,2,2) 
scatter3(p_L2(3,:),p_L2(2,:),p_L2(1,:),2,sqrt(sum(p_S2.^2,1))) 
axis vis3d 


nan(0); 

