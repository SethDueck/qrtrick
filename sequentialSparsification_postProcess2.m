

filterb3 = false([size(ssdata,1) 1]);
filterb3( indexListRound2(filter2) ) = true; 
filterb3(filter1) = true; 
indexListRound3 = find(~filterb3); 
ssdata3 = ssdata(~filterb3,:); 

%
% Mostly copied from "first run" code 
%
% ssdata2 = ssdata2 -  repmat(mean(ssdata2,  1), [size(ssdata2,1),1]); 
% ssdata2 = ssdata2 ./ repmat(std( ssdata2,0,1), [size(ssdata2,1),1]); 

[mu3,lambda3] = getDefaultRpcaParams(ssdata3); 

figure(15) 
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
filter3 = 1:size(ssdata3,1); 
filter3_smallSparsePass = true([size(ssdata3,1) 1]); 
S3_old = nan(1);  
Y3_old = nan(1); 
isConverged = false; 
LRank3_old = []; 
for erroridx = 1:numel(errorSequence)
    if( isConverged ) 
        break; 
    end 
    
    filter3 = filter3(filter3_smallSparsePass); 
%     [L1,S1] = rpca(ssdata(filter,:),1e-2, mu, lambda);
    [mu3,lambda3] = getDefaultRpcaParams(ssdata3(filter3,:)); 
    [L3,S3, Y3, LRank3] = rpca(ssdata3(filter3,:),errorSequence(erroridx), mu3, lambda3, S3_old, Y3_old, LRank3_old ); 
    % X = L1 + S1 + Y/mu -- this *may* be true after the first pass 

    unsmoothPart3 = ssdata3(filter3,:) - L3; 
    unsmoothPart3_mag = sqrt(sum(unsmoothPart3.^2,2)); 
    filter3_smallSparsePass = (unsmoothPart3_mag-min(unsmoothPart3_mag)) < ... 
        ((max(unsmoothPart3_mag)-min(unsmoothPart3_mag))/4); 
    
    fracErrorL  = sqrt(sum((ssdata3(filter3,:)-L3   ).^2,2)) ...
        ./ sqrt(sum((ssdata3(filter3,:)).^2,2)); 
    fracErrorLS = sqrt(sum((ssdata3(filter3,:)-L3-S3).^2,2)) ...
        ./ sqrt(sum((ssdata3(filter3,:)).^2,2));     
    isConverged = max(fracErrorLS) < 1e-3; 
    
    
	Y3_old = Y3(filter3_smallSparsePass,:); 
    S3_old = S3(filter3_smallSparsePass,:); 
    L3_old = L3; 
    LRank3_old = LRank3; 
    
    figure(15) 
    subplot(4,1,1)
    plot(filter3, sqrt(sum(S3.^2,2)) ,'.')
    title([num2str(erroridx) ': ' num2str(numel(filter3))]); 
    ylabel('S') 
    subplot(4,1,2) 
    plot(filter3, sqrt(sum((ssdata3(filter3,:)-L3).^2,2)) ,'.')
    ylabel('data - L') 
    nan(0); 
    subplot(4,1,3) 
    plot(filter3, sqrt(sum((ssdata3(filter3,:)-L3-S3).^2,2)) ,'.')
    ylabel('data - L - S') 
    subplot(4,1,4)
    cla 
    hold all 
    plot(filter3, fracErrorL ,'.')
    plot(filter3, fracErrorLS,'.') 
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


[u_L3,s_L3,v_L3] = svd(L3,'econ'); 
[u_S3,s_S3,v_S3] = svd(S3,'econ'); 
[u_r3,s_r3,v_r3] = svd(ssdata3(filter3,:)-L3,'econ'); 

p_L3 = v_L3'*L3'; 
p_S3 = v_S3'*S3'; 
p_r3 = v_r3'*(ssdata3(filter3,:)-L3)'; 

hr_subset = hr(~filterb3,:); 
figure(13) 
clf 
subplot(1,2,1) 
scatter3(p_L3(3,:),p_L3(2,:),p_L3(1,:),1,hr_subset(filter3,10))
axis vis3d 
subplot(1,2,2) 
scatter3(p_L3(3,:),p_L3(2,:),p_L3(1,:),2,sqrt(sum(p_S3.^2,1))) 
axis vis3d 


nan(0); 

