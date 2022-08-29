
%
% This is the original "sandbox" version of my sequential sparsification
% algorithm. It was turned into a function "sequentialSparsification_step"
% that is run iteratively by a script "sequentialSparsification" that
% handles loading the data and passing subsets of data to calls of the
% function. 
% 


% Use a "typical" normalized ssalbedo dataset for prototyping 
% ssdata = data; 
if(~exist('ssdata','var')) 
    temp = load('ssdataset.mat'); 
    ssdata = temp.data; 
%     warning('Haven''t acutally loaded this dataset before, probably have to pull data out of a temp var.'); 
    
end 

ssdata = ssdata -  repmat(mean(ssdata,  1), [size(ssdata,1),1]); 
ssdata = ssdata ./ repmat(std( ssdata,0,1), [size(ssdata,1),1]); 

[mu,lambda] = getDefaultRpcaParams(ssdata); 

figure(4) 
clf 
subplot(4,1,1)
hold all 
subplot(4,1,2)
hold all 
subplot(4,1,3)
hold all 

figure(6) 
clf 

errorSequence = 10.^linspace(-2,-5,20); 
% errorSequence = 10.^linspace(-2,-2,30); 
% filter = true([size(ssdata,1) 1]); 
filter1 = 1:size(ssdata,1); 
filter_smallSparsePass = true([size(ssdata,1) 1]); 
isConverged = false; 
S1_old = nan(1);  
Y1_old = nan(1); 
tic
timevars = nan([numel(errorSequence),8]); 
LRank1_old = []; 
for erroridx = 1:numel(errorSequence)
    if( isConverged ) 
        break; 
    end 
    
    if(10==erroridx) 
       nan(0); 
    end 

    filter1 = filter1(filter_smallSparsePass); 
    [mu,lambda] = getDefaultRpcaParams(ssdata(filter1,:)); 
    convergenceFactor = errorSequence(erroridx); 
%     if( 1==erroridx ) 
%         convergenceFactor = 1e-2; 
%     else 
%         convergenceFactor = min(fracErrorL); 
%     end

tic; 
    [L1,S1, Y1, LRank1] = rpca(ssdata(filter1,:),convergenceFactor, mu, lambda, S1_old, Y1_old, LRank2_old ); 
etime = toc;     
    
    unsmoothPart = ssdata(filter1,:) - L1; 
    unsmoothPart_mag = sqrt(sum(unsmoothPart.^2,2)); 
    filter_smallSparsePass = (unsmoothPart_mag-min(unsmoothPart_mag)) < ... 
        ((max(unsmoothPart_mag)-min(unsmoothPart_mag))/4); 
    
    fracErrorL  = sqrt(sum((ssdata(filter1,:)-L1   ).^2,2)) ...
        ./ sqrt(sum((ssdata(filter1,:)).^2,2)); 
    fracErrorLS = sqrt(sum((ssdata(filter1,:)-L1-S1).^2,2)) ...
        ./ sqrt(sum((ssdata(filter1,:)).^2,2));     
    isConverged = max(fracErrorLS) < 1e-3; 
    
%     figure(6) 
%     if(1~=erroridx) 
%     subplot(2,6,erroridx) 
%     hold all 
% %     plot(sort(S1_old(:)-S1(:))); 
% %     temp1 = S1_old(filter_smallSparsePass,:);
% %     temp2 = S1(    filter_smallSparsePass,:); 
% %     plot(sort(temp1(:)-temp2(:))); 
%     [tempy,tempx] = hist(S1_old(:)-S1(:),1001); 
%     tempy(tempy==max(tempy)) = 0; 
%     plot(tempx,tempy) 
%     title([num2str(sum(abs(S1_old(:)-S1(:))>1e-2) ) ': ' num2str(etime)]) 
%     end
    
    Y1_old = Y1(filter_smallSparsePass,:); 
    S1_old = S1(filter_smallSparsePass,:); 
    L1_old = L1; 
    LRank1_old = LRank1; 
    
    figure(4) 
    subplot(4,1,1)
    plot(filter1, sqrt(sum(S1.^2,2)) ,'.')
    title([num2str(erroridx) ': ' num2str(numel(filter1))]); 
    ylabel('S') 
    subplot(4,1,2) 
    plot(filter1, sqrt(sum((ssdata(filter1,:)-L1).^2,2)) ,'.')
    ylabel('data - L') 
    nan(0); 
    subplot(4,1,3) 
    plot(filter1, sqrt(sum((ssdata(filter1,:)-L1-S1).^2,2)) ,'.')
    ylabel('data - L - S') 
    subplot(4,1,4)
    cla 
    hold all 
    plot(filter1, fracErrorL ,'.')
    plot(filter1, fracErrorLS,'.') 
    ylabel('frac error') 
    title(['maxL = ' num2str(max(fracErrorL )) ',   maxLS = ' ... 
        num2str(max(fracErrorLS)) ',   e = ' num2str(convergenceFactor)]); 
    
    
    timevars(erroridx,1) = etime; 
    timevars(erroridx,2) = convergenceFactor; 
    timevars(erroridx,3) = min(fracErrorL); 
    timevars(erroridx,4) = max(fracErrorL); 
    timevars(erroridx,5) = min(fracErrorLS); 
    timevars(erroridx,6) = max(fracErrorLS); 
    timevars(erroridx,7) = mu; 
    timevars(erroridx,8) = lambda; 
%     figure(5) 
%     clf 
%     plot3(sqrt(sum((ssdata(filter,:)-L1-S1).^2,2)), sqrt(sum((ssdata(filter,:)-L1).^2,2)), sqrt(sum(S1.^2,2)),'.') 
%     xlabel('data - L - S') 
%     ylabel('data - L') 
%     zlabel('S') 
%     axis vis3d

    
    nan(0); 
end
toc 

[u_L1,s_L1,v_L1] = svd(L1,'econ'); 
[u_S1,s_S1,v_S1] = svd(S1,'econ'); 
[u_r1,s_r1,v_r1] = svd(ssdata(filter1,:)-L1,'econ'); 

p_L1 = v_L1'*L1'; 
p_S1 = v_S1'*S1'; 
p_r1 = v_r1'*(ssdata(filter1,:)-L1)'; 

figure(2) 
clf 
subplot(1,2,1) 
scatter3(p_L1(3,:),p_L1(2,:),p_L1(1,:),2,hr(filter1,40))
axis vis3d 
subplot(1,2,2) 
scatter3(p_L1(3,:),p_L1(2,:),p_L1(1,:),2,sqrt(sum(p_S1.^2,1)))
axis vis3d 




nan(0); 

