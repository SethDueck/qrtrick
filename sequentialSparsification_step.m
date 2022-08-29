function [ filter_wasClustered, L, LRank, S, timedata ] = sequentialSparsification_step( ssdata, errorSequence, maxLRank, convergenceLSError, srdopts )


if( ~isnan(srdopts.fidx_iterationConvergence) )
    figure(srdopts.fidx_iterationConvergence) 
    clf 
%     subplot(4,1,1)
%     hold all 
%     subplot(4,1,2)
%     hold all 
%     subplot(4,1,3)
%     hold all 
end 

% Set initial guesses to nan for variables that get used between iterations 
filter_localIndices_old = 1:size(ssdata,1); 
filter_smallSparsePass = true([size(ssdata,1) 1]); 
S_old = nan(1);  
Y_old = nan(1); 
LRank_old = []; 


% [mu,lambda] = getDefaultRpcaParams(ssdata(filter_localIndices,:)); 

    
timedata = nan([numel(errorSequence),8]); 
isConverged = false; 
for erroridx = 1:numel(errorSequence)
    if( isConverged ) 
        break; 
    end 
    
    if(10==erroridx) 
       nan(0); 
    end 
    
    
    % Do robust PCA (sparsification) step 
    % Skip step if data is already low-rank 
    % Get "a priori" (before sparsification) eigenvalues 
%     [vAP,dAP] = eig(ssdata(filter_localIndices_old,:)'*ssdata(filter_localIndices_old,:)); 
    [ut,dAP,vt] = srd_svd(ssdata(filter_localIndices_old,:),0); 
    thresh_sparsificationNumDimensions = 5; 
    isNeedSparsification = dAP(thresh_sparsificationNumDimensions)>1; 
    
    for uidx = 1:5 
        % Roughly calculate sparsest part of distribution of ut(:,uidx); 
        % take midpoint of sparsest area as split point
        [un,ux] = hist(ut(:,uidx)); 
        [~,midx] = min(un); 
        utSplitPoint = ux(midx); 
%         uts = sort(ut(:,uidx)); 
%         [~,midx] = max(diff(uts)); 
%         utSplitPoint = 0.5*(uts(midx)+uts(midx+1)); 
        std_whole = std(ut(:,uidx)); 
        std_pos = std(ut(ut(:,uidx) >= utSplitPoint,uidx)); 
        std_neg = std(ut(ut(:,uidx) <  utSplitPoint,uidx)); 
        if( std_whole > 2*sqrt(std_pos*std_pos+std_neg*std_neg) ) 
            % Split on ut 
            % Keep larger subset 
            if(sum(ut(:,uidx)>utSplitPoint) > 0.5*size(ut,1)) 
                filter_ut = ut(:,uidx)>utSplitPoint; 
            else 
                filter_ut = ut(:,uidx)<utSplitPoint; 
            end 
            L_old = L(filter_ut,:); 
            S_old = S(filter_ut,:); 
            Y_old = Y(filter_ut,:); 
            filter_localIndices_old = filter_localIndices_old(filter_ut); 
            break; % Only filter on one ut dimension -- just to avoid writing logic for splitting on multiple dimensions 
        end 
    end 
    
    convergenceFactor = errorSequence(erroridx); 
        
    if( isNeedSparsification )
    [eigShrinkMagnitude,lambda] = getDefaultRpcaParams(ssdata(filter_localIndices_old,:)); 

    % Split data into Low-rank + Sparse matrices 
    tic; 
    if(LRank_old>maxLRank) 
        maxNextLRank = LRank_old; 
    else 
        maxNextLRank = LRank_old; 
    end 
%     maxNextLRank = 0; 
    [L,S, Y, LRank, U_last, S_reduced] = rpca(ssdata(filter_localIndices_old,:),convergenceFactor, eigShrinkMagnitude, lambda, S_old, Y_old, maxNextLRank );
%     [L,S, Y, LRank] = rpca(ssdata(filter_localIndices_old,:),convergenceFactor, mu, lambda, [], [], LRank_old ); 
    etime = toc;     
    
    % Find where PCs contribute sparsely to reconstruction 
    numUBetweenKnees = nan([LRank 1]); 
    numUOutsideKnees = nan([LRank 2]); 
    sqSumOutsideKnees = nan([LRank 2]); 
    stdKnees = nan([LRank 3]); 
    meanKnees = nan([LRank 3]); 
    kneeSlopes = nan([LRank 3]); 
    uKnees = zeros([size(U_last,1) LRank]); 
    filter_slopedTail = true([size(U_last,1) 1]); 
    U_last_sig = U_last(:,1:LRank); 
    for uidx = 1:LRank
        numUBetweenKnees(uidx) = findNumBetweenKnees(U_last(:,uidx));
        f_p = (U_last(:,uidx) >=  findKnee( U_last(:,uidx))); 
        f_n = (U_last(:,uidx) <= -findKnee(-U_last(:,uidx))); 
        numUOutsideKnees(uidx,2) = sum(f_p); 
        sqSumOutsideKnees(uidx,2) = sum(U_last(f_p,uidx).^2); 
        
        numUOutsideKnees(uidx,1) = sum(f_n); 
        sqSumOutsideKnees(uidx,1) = sum(U_last(f_n,uidx).^2); 
        
        stdKnees(uidx,1) = std(U_last(f_n,uidx)); 
        stdKnees(uidx,2) = std(U_last(~(f_n|f_p),uidx)); 
        stdKnees(uidx,3) = std(U_last(f_p,uidx)); 
        
%         ytotal = max(U_last(:,uidx))-min(U_last(:,uidx)); 
%         ytotal = (max(U_last(~(f_n|f_p),uidx))-min(U_last(~(f_n|f_p),uidx))); 
%         slopex = numUOutsideKnees(uidx,1)/size(U_last,1); 
%         slopey = (max(U_last(f_n,uidx))-min(U_last(f_n,uidx)))/ytotal; 
%         kneeSlopes(uidx,1) = slopey/slopex; 
%         slopex = numUBetweenKnees(uidx)/size(U_last,1); 
%         slopey = (max(U_last(~(f_n|f_p),uidx))-min(U_last(~(f_n|f_p),uidx)))/ytotal; 
%         kneeSlopes(uidx,2) = slopey/slopex; 
%         slopex = numUOutsideKnees(uidx,2)/size(U_last,1); 
%         slopey = (max(U_last(f_p,uidx))-min(U_last(f_p,uidx)))/ytotal; 
%         kneeSlopes(uidx,3) = slopey/slopex; 
%         uKnees(f_n|f_p,uidx) = U_last(f_n|f_p,uidx); 

        [ulim_n,ulim_p,kneeBadnessMeasure] = findKnee_subSlope(U_last(:,uidx)); 
        kneeSlopes(uidx,:) = kneeBadnessMeasure; 
        if((kneeBadnessMeasure(2)<kneeBadnessMeasure(1)) && (kneeBadnessMeasure(2)<kneeBadnessMeasure(3)) )
            filter_temp = (U_last_sig(:,uidx)>ulim_n) & (U_last_sig(:,uidx)<ulim_p); 
            U_last_sig(filter_temp,uidx) = 0; 
        else 
            U_last_sig(:,uidx) = 0; 
        end
        
        figure(311)
        clf
        hold all
        plot(sort(U_last(:,uidx)))
        plot([1 size(U_last,1)],[1 1]*ulim_n)
        plot([1 size(U_last,1)],[1 1]*ulim_p)

%         kneeSlopeLimit = 20; 
%         if(kneeSlopes(uidx,1)>kneeSlopeLimit) 
%             filter_slopedTail = filter_slopedTail & ~f_n; 
%         end
%         if(kneeSlopes(uidx,3)>kneeSlopeLimit) 
%             filter_slopedTail = filter_slopedTail & ~f_p; 
%         end 
        
    end
    U_last_sig = U_last_sig*S_reduced(1:LRank,1:LRank); 
    [vu,du] = eig(U_last_sig'*U_last_sig);
pu = vu'*U_last_sig';
pu_mag = sqrt(sum(pu.^2,1));
% filter_slopedTail = pu_mag==0; 
filter_slopedTail = pu_mag < findKnee(pu_mag); 
    filter_localIndices_old = filter_localIndices_old(filter_slopedTail); 
    L = L(filter_slopedTail,:); 
    S = S(filter_slopedTail,:); 
    Y = Y(filter_slopedTail,:); 

    
    % Filter out wavelengths with a large sparse component 
    S_mag = sqrt(sum(S.^2,2)); 
    filter_S = S_mag <= findKnee(S_mag); 
    filter_localIndices_old = filter_localIndices_old(filter_S); 
    L = L(filter_S,:); 
    S = S(filter_S,:); 
    Y = Y(filter_S,:); 
    else 
        L = L_old; 
        S = S_old; 
        
        % Need to update L, S, to reflect removal of wavelengths filtered
        % out in last round 
        [vd,dd] = eig(ssdata(filter_localIndices_old,:)'*ssdata(filter_localIndices_old,:));
        [~,sortfilter] = sort(real(diag(dd)), 'descend'); 
        vd = vd(:,sortfilter); 
        vL = vd(:,1:thresh_sparsificationNumDimensions); 
        L = (ssdata(filter_localIndices_old,:)*vL)*vL';
        S = shrink(ssdata(filter_localIndices_old,:)-L, lambda*eigShrinkMagnitude);
    end 
    
    
    % Find wavelengths that are least-well-described by current low-rank
    % factorization 
    unsmoothPart = ssdata(filter_localIndices_old,:) - L-S; 
    unsmoothPart_mag = sqrt(sum(unsmoothPart.^2,2)); 
%     unsmoothPart_limit = (max(unsmoothPart_mag)-min(unsmoothPart_mag)) / 4; 
    unsmoothPart_limit = findKnee(unsmoothPart_mag); 
    filter_smallSparsePass = (unsmoothPart_mag-min(unsmoothPart_mag)) < ...
        unsmoothPart_limit; 
    
    % Identify principal components "specific" to those wavelengths 
    [v1,s1] = eig(unsmoothPart'*unsmoothPart); 
    [s1,sortfilter] = sort(real(diag(s1)),1,'descend'); 
    v1 = v1(:,sortfilter); 
%     [u2,s2,v2] = svd(unsmoothPart( filter_smallSparsePass,:),'econ'); 
%     [u3,s3,v3] = svd(,'econ'); 
    [v3,s3] = eig(unsmoothPart(~filter_smallSparsePass,:)'*unsmoothPart(~filter_smallSparsePass,:)); 
    [s3,sortfilter] = sort(real(diag(s3)),1,'descend'); 
    v3 = v3(:,sortfilter); 
    isUnsmoothCompImportant = (s3/sum(s3) > s1/sum(s1)) | ((s3/sum(s3))>0.1); 
    isUnsmoothCompImportant = isUnsmoothCompImportant & (s3/sum(s3) > 1e-2); 
    isUnsmoothCompImportant(1) = true; 
%     s1_iuci = s1; 
%     s1_iuci(s1_iuci<0) = 0; 
%     s3_iuci = s3; 
%     s3_iuci(s3_iuci<0) = 0; 
%     isUnsmoothCompImportant = ( (s3_iuci/sum(s3_iuci)) > (s1_iuci/sum(s1_iuci)) ) | ... 
%         ((s3_iuci/sum(s3_iuci))>0.1); 
%     isUnsmoothCompImportant = isUnsmoothCompImportant & ((s3_iuci/sum(s3_iuci))>1e-5); 
%     isUnsmoothCompImportant(end) = false; 
%     tempidx = find(isUnsmoothCompImportant,1,'last'); 
%     if(10<tempidx)
%         warning('large number of unsmooth components to be removed'); 
%     end 
%     isUnsmoothCompImportant(1:tempidx) = true; 
    vUnsmooth = v3(:,isUnsmoothCompImportant); 
    
    % Identify wavelengths that have those components 
    pUnsmooth = vUnsmooth'*unsmoothPart'; 
    pUnsmooth_mag = sqrt(sum(pUnsmooth.^2,1)); 
    pUnsmooth_limit = findKnee(pUnsmooth_mag);     
    filter_pUnsmooth = pUnsmooth_mag < pUnsmooth_limit; 
    filter_pUnsmoothComponents = true(size(filter_pUnsmooth)); 
%     for pucidx = 1:size(pUnsmooth,1) 
%         filter_pUnsmoothComponents = filter_pUnsmoothComponents & ... 
%           (  pUnsmooth(pucidx,:) < findKnee( pUnsmooth(pucidx,:)) ) & ...
%           ( -pUnsmooth(pucidx,:) < findKnee(-pUnsmooth(pucidx,:)) ); 
%     end 
    filter_pUnsmooth = filter_pUnsmooth & filter_pUnsmoothComponents; 
    filter_smallSparsePass = filter_pUnsmooth'; 
    
    % TESTING: Don't check for unsmooth components 
    filter_smallSparsePass = true(size(filter_smallSparsePass)); 
    
    % Remove those wavelengths from this cluster 
    filter_localIndices = filter_localIndices_old(filter_smallSparsePass); 
    S = S- (S*vUnsmooth)*vUnsmooth'; 
    Y = Y- (Y*vUnsmooth)*vUnsmooth'; 
    
    % Calculate error between "target" and L+S decomposition 
    fracErrorL  = sqrt(sum((ssdata(filter_localIndices_old,:)-L   ).^2,2)) ...
        ./ sqrt(sum((ssdata(filter_localIndices_old,:)).^2,2)); 
    fracErrorLS = sqrt(sum((ssdata(filter_localIndices_old,:)-L-S).^2,2)) ...
        ./ sqrt(sum((ssdata(filter_localIndices_old,:)).^2,2));     
    isConverged = (maxLRank >= LRank) && ... 
        max(fracErrorLS(filter_smallSparsePass)) < convergenceLSError; 
    
    % Update speedup variables for wavelengths that pass this iteration 
    L_old = L(filter_smallSparsePass,:); 
    S_old = S(filter_smallSparsePass,:); 
    Y_old = Y(filter_smallSparsePass,:); 
    LRank_old = LRank; 
    
    if( ~isnan(srdopts.fidx_iterationConvergence) )
        
        figure(srdopts.fidx_iterationConvergence+1) 
        [~,~,vL] = srd_svd(L,0); 
        plot(vL(:,1:3)); 
        
        figure(srdopts.fidx_iterationConvergence) 
        subplot(4,1,1)
        cla 
        hold all 
        plot(filter_localIndices_old, sqrt(sum(S.^2,2)) ,'.')
        plot(filter_localIndices, sqrt(sum(S(filter_smallSparsePass,:).^2,2)),'.')
        title([num2str(erroridx) ': ' num2str(sum(filter_smallSparsePass)) ' / ' num2str(numel(filter_localIndices_old)) '   (' num2str(size(ssdata,1)) '),   rank=' num2str(LRank)]); 
        ylabel('S') 
        subplot(4,1,2) 
        cla
        hold all 
        plot(filter_localIndices_old, unsmoothPart_mag ,'.')
        plot(filter_localIndices,unsmoothPart_mag(filter_smallSparsePass),'.')
        ylabel('data - L') 
        title(['knee L: ' num2str(findKnee(fracErrorL))]); 
        nan(0); 
        subplot(4,1,3) 
        cla 
        hold all 
        plot(filter_localIndices_old, sqrt(sum((ssdata(filter_localIndices_old,:) - L - S ).^2,2)) ,'.')
        plot(filter_localIndices, sqrt(sum((ssdata(filter_localIndices,:)-L(filter_smallSparsePass,:)-S(filter_smallSparsePass,:)).^2,2)),'.')
        ylabel('data - L - S') 
        title(['knee LS: ' num2str(findKnee(fracErrorLS))]); 
        subplot(4,1,4)
        cla 
        hold all 
        plot(filter_localIndices_old, fracErrorL ,'.')
        plot(filter_localIndices_old, fracErrorLS,'.') 
        ylabel('frac error') 
        title(['maxL = ' num2str(max(fracErrorL )) ',   maxLS = ' ... 
            num2str(max(fracErrorLS)) ',   e = ' num2str(convergenceFactor)]); 
    end
    
    filter_localIndices_old = filter_localIndices; 
    
    % Keep track of simple timing statistics 
    timedata(erroridx,1) = etime; 
    timedata(erroridx,2) = convergenceFactor; 
    timedata(erroridx,3) = min(fracErrorL); 
    timedata(erroridx,4) = max(fracErrorL); 
    timedata(erroridx,5) = min(fracErrorLS); 
    timedata(erroridx,6) = max(fracErrorLS); 
    timedata(erroridx,7) = eigShrinkMagnitude; 
    timedata(erroridx,8) = lambda; 
    
    nan(0); 
end

L = L(filter_smallSparsePass,:); 
S = S(filter_smallSparsePass,:); 
filter_wasClustered = filter_localIndices; 


nan(0); 


end 
