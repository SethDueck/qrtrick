

numQRVecs = 20; 
qrdata = ssdata; 

qrvecs   = nan(size(qrdata,1),numQRVecs); 
qrproj = nan(numQRVecs,size(qrdata,2)); 

for nidx = 1:numQRVecs 
    [~,maxidx] = max(sum(abs(qrdata),1)); 
    
    qrunit = qrdata(:,maxidx); 
    qrunit = qrunit / sqrt(sum(qrunit.^2)); 
%     qrscores = qrdata*qrunit'; 
    qrscores = qrunit'*qrdata; 
    
    qrvecs(:,nidx) = qrunit; 
    qrproj(nidx,:) = qrscores; 
    qrdata = qrdata - qrunit*qrscores; 
    
    nan(0); 
    
end 

qrproj_norm = qrproj; 
for didx = 1:size(qrproj_norm,1) 
    qrproj_norm(didx,:) = qrproj_norm(didx,:) / sqrt(sum(qrproj_norm(didx,:).^2)); 
end 

[v2,d2] = eig(abs(qrvecs)'*abs(qrvecs) - qrvecs'*qrvecs); 
[d2,sortfilter] = sort(diag(real(d2)),'descend'); 
v2 = v2(:,sortfilter); 
a = v2'*qrproj; 

% Find most linearly dependent variables 
a_n = a; 
for aidx = 1:size(a,1) 
    a_n(aidx,:) = a_n(aidx,:) / sqrt(sum(a_n(aidx,:).^2)); 
end 
dataproj = ssdata*a_n'; 
% scatter3(temp(:,end-2),temp(:,end-1),temp(:,end),2,fs(:,70)./ss(:,70))
dCovMatrix = nan([1 1]*size(dataproj,2)); 
for bidx1 = 1:size(dCovMatrix,1) 
    dCovMatrix(bidx1,bidx1) = 0; % Give value to diagonal so this can be plotted 
    for bidx2 = (bidx1+1):size(dCovMatrix,2) 
        dCovMatrix(bidx1,bidx2) = det(dataproj(:,[bidx1 bidx2])'*dataproj(:,[bidx1 bidx2])); 
        dCovMatrix(bidx2,bidx1) = dCovMatrix(bidx1,bidx2); 
    end 
end 

% Sort other dimensions bidx2 according to how linear dependent they are
% with bidx1
bmax = max(max(abs(dCovMatrix))); 
mostLinearDep = nan([1 1]*size(dCovMatrix,2)); % Each row is sorted, best candidate for xidx is at yidx=1 
linearDepMeasure = nan([1 1]*size(dCovMatrix,2)); 
for bidx = 1:size(dCovMatrix,2) 
    dCovMatrix(bidx,bidx) = 1.1*bmax; % Set diagonal to large value so we can find minimum of non-zero values 
    [linearDepMeasure(bidx,:),mld] = sort(dCovMatrix(bidx,:)); 
    mostLinearDep(bidx,:) = mld; 
end 
% plot(sum(linearDepMeasure(:,1:4),2))
[~,linearDepOrder] = sort(linearDepMeasure(:,1)); 

% Use sindy to find order 1 relationships 
pointWeights = ones([size(dataproj,1) 1]); 
maxNumLinTerms = 6; 
for linidx = 1:maxNumLinTerms  
    if(linidx > size(dataproj,2)) 
        warning('Trying to fit to more dims than there are dims remaining'); 
        break; 
    end 
    
    if(numel(linearDepOrder)~=size(dataproj,2)) 
        % Lazy way to avoid recalculating linearDepOrder after the first
        % round. Not sure if this the order is actually useful anyway;
        % either replace this with real code or delete linearDepOrder. 
        linearDepOrder = (1:size(dataproj,2))'; 
    end 
    order1Error = nan(size(linearDepOrder)); 
    regressionModel = nan([1 1]*size(dataproj,2) + [0 1]); 
    
for vidx = 1:size(dataproj,2) 
    
    dimidx = linearDepOrder(vidx); 
    
%     A = vmean; 
%     b = dvdt(:,dimidx); 
    if( (dimidx>size(mostLinearDep,1)) || (linidx>size(mostLinearDep,2)) ) 
        nan(0); 
    end 
    indepDimIdx = mostLinearDep(dimidx,1:linidx); 
    A = [ dataproj(:, indepDimIdx), ... % Treat most linear dimensions as independent 
          ones([size(dataproj,1) 1]); ]; % Add in affine transformation even though the zeroth order term seems to always be zero atm 
    b = dataproj(:,dimidx); 
    C = eye(size(A,2)); 
    
    eta = 1e2; 
    lambda = 1e-1; 
    w0 = ones([size(A,2) 1]); 
    convergenceThreshold = 1e-5; 
    maxNumIterations = 1e4; 
    [ xk, wk ] = srdsr3_mod( A, C, b, pointWeights, eta, lambda, w0, convergenceThreshold, maxNumIterations ); 
    
    % Check fit 
    order1Fit = A*xk; 
    order1Error(dimidx) = sqrt(sum((b-order1Fit).^2)); 
    
    regressionModel(dimidx,:) = 0; 
    regressionModel(dimidx,indepDimIdx) = xk(1:numel(indepDimIdx)); 
    regressionModel(dimidx,dimidx) = -1; 
    regressionModel(dimidx,end) = xk(end); 
    
    nan(0); 
    
end 

fracError = order1Error ./ sqrt(sum(dataproj.^2,1)'); 

% Normalize equations 
normalizedModel = regressionModel;
for nidx = 1:size(regressionModel,1)
    normalizedModel(nidx,:) = normalizedModel(nidx,:) / sqrt(sum(normalizedModel(nidx,:).^2));
end

% Look for degeneracy in system 
modelDeterminate = abs(det(regressionModel(:,1:end-1))); 
isDegenerate =  modelDeterminate < 1e-5; 

% Take degeneracy to mean that algorithm has converged on equations --
% reduce number of rows in system by replacing degenerate equations with an
% auxiliary coordinate 
% BUT fit must also be good (check that order1Error is low) 

% Eigenvectors with eigenval>1 tell which variables are dependent: 
modelCov = abs(normalizedModel*normalizedModel'); 
modelCov_bound = 0.99; 
maxOffDiagonal = max(reshape(triu(modelCov,1)+tril(modelCov,-1),1,[])); 
if(modelCov_bound > maxOffDiagonal) % Get at least one pair of degenerate equations 
    modelCov_bound = maxOffDiagonal; 
end
modelCov_thresholded = (abs(normalizedModel*normalizedModel')>=modelCov_bound)+0; 
[v_ideal,d_ideal] = eig(modelCov_thresholded); 
[d_ideal, sortfilter] = sort(real(diag(d_ideal)),'descend'); 
v_ideal = v_ideal(:,sortfilter); 

% Normalize projected data matrix 
dpn = dataproj;
for tempidx = 1:size(dpn,2) 
dpn(:,tempidx) = dpn(:,tempidx) / sqrt(sum(dpn(:,tempidx).^2));
end

% Find degenerate dimensions from set of equations 
numDegenSets = sum((d_ideal>1)+0); 
degenColumns = cell([numDegenSets 1]); 
degenError   = nan([ numDegenSets 1]); 
for didx = 1:numDegenSets; 
    degenerateEquations = find(abs(v_ideal(:,didx)) > 1e-5); % These equations are degenerate, 
    thisSetDegenCols = find( sum(abs(regressionModel(degenerateEquations,:)),1) > 1e-5); % so collect a linear combination of terms from them. 
    degenColumns{didx} = thisSetDegenCols; 
    degenError(didx)   = max(fracError(thisSetDegenCols)); 
end
[~,degenSortOrder] = sort(degenError); % Order degeneracies best to worst fit 

% Consolidate dependent variables 
columnsToDelete = []; 
touchedColumns = []; 
for didx = 1:numel(degenSortOrder) 

    % Recall degenerate columns 
    thisSetDegenCols = degenColumns{degenSortOrder(didx)}; 
    degdet = det(dataproj(:,thisSetDegenCols)'*dataproj(:,thisSetDegenCols)); 
    theseEqnErrors = order1Error(thisSetDegenCols); % Did best model actually fit the data? 
    
    % Save relationship between variables and new auxiliary variable xi
    % TODO 
    
    % Project onto pcs of degenerate space 
    [v_deg,d_deg] = eig(dataproj(:,thisSetDegenCols)'*dataproj(:,thisSetDegenCols)); 
    [d_deg,sortfilter] = sort(real(diag(d_deg)),'descend'); 
    v_deg = v_deg(:,sortfilter); 
    p_deg = v_deg'*dataproj(:,thisSetDegenCols)'; 
    % OR do QR pivot trick to generate "mod-PCs" 
    degenData = dataproj(:,thisSetDegenCols); 
    degenData_qr = nan(size(degenData)); 
    for ddidx = 1:size(degenData,2) 
        [~,maxidx] = max(sum(abs(degenData),2)); 
        unitv = degenData(maxidx,:); 
        unitv = unitv / sqrt(sum(unitv.^2)); 
        unitvscores = degenData*unitv'; 
        degenData_qr(:,ddidx) = unitvscores; 
        degenData = degenData - unitvscores*unitv; 
    end 
    qrInfoShare = sum(degenData_qr'.^2,2)/sum(sum(degenData_qr.^2)); 
    
    % Decide which dims to keep, replace data in dataproj 
    keepDimFlag = (d_deg/sum(d_deg)) > 0.01; 
%     keepDimFlag(end) = false; % Force deletion of least important dim 
    if( ~keepDimFlag(end) ) 
        % Replace dims with pcs if data loss is within tolerance 
        touchedColumns  = [touchedColumns;  thisSetDegenCols( keepDimFlag)'];  %#ok<AGROW>
        columnsToDelete = [columnsToDelete; thisSetDegenCols(~keepDimFlag)'];  %#ok<AGROW>
        dataproj(:,thisSetDegenCols(keepDimFlag)) = p_deg(keepDimFlag,:)'; 
%         replacementColumns = [replacementColumns; ones([sum(~keepDimFlag+0) 1])*thisSetDegenCols(find(keepDimFlag,1,'first'))];  %#ok<AGROW>
    end 
    nan(0); 

    % Replace variables with xi 
    % For now I will just exclude columns -- should really create var xi
    % and replace ALL data with "smoothed" data, not just remove one
    % dimension at random. 
    % AND the actual deletion happens later -- I'm really just making a
    % list of what to delete later. 
    % Decide how many columns to exclude 

    nan(0); 
    
end 

if( any( ismember(columnsToDelete,touchedColumns) ) ) 
    warning('Replacements and to-delete overlaps') 
end 
if( numel(columnsToDelete) ~= numel(unique(columnsToDelete)) ) 
    warning('Columns deleted more than once') 
end 

% % Delete columns 
% % Change references to columns that will be deleted 
% for cidx = 1:numel(columnsToDelete) 
%     mostLinearDep(mostLinearDep(:)==columnsToDelete(cidx)) = replacementColumns(cidx); % Update indexing for references to rows that will be deleted. Should be updated to reference xi; currently I just ref back to the not-deleted data column 
% end 
% % Check for self reference -- we will get x_i=x_i from the regression if
% % x_i is fed in as an independent variable on the next round. 
% mayHaveSelfReference = true; 
% while( mayHaveSelfReference ) 
%     mayHaveSelfReference = false; % Repeat in case self-reference is found 
% for cidx = 1:size(dataproj,2) 
%     if( cidx==mostLinearDep(cidx,1) ) 
%         % Bump self-reference to the end of the list of candidate indep 
%         % vars; just permute the sequence. There may be problems because I
%         % don't check for uniqueness within the rows here. 
%         mostLinearDep(cidx,1:end-1) = mostLinearDep(cidx,2:end); 
%         mostLinearDep(cidx,end)     = mostLinearDep(cidx,end); 
%         mayHaveSelfReference = true; 
%         nan(0); 
%     end 
% end 
% end 
% 
% % Indexing will change when dims are deleted -- figure out how much to subtract from each dimension index 
% dimIndexShift = cumsum( ismember(1:size(dataproj,2),columnsToDelete) + 0 ); 
% % Do the subtraction 
% for cidx = 1:size(dataproj,2) 
%     mostLinearDep(mostLinearDep(:)==cidx) = cidx - dimIndexShift(cidx); 
% end 
% mostLinearDep(columnsToDelete,:) = []; % Delete rows because dimensions are being deleted 
% % Delete redundant columns 
% mostLinearDep_new = nan([1 1]*size(mostLinearDep,1)); 
% for cidx = 1:size(mostLinearDep,1) 
%     mostLinearDep_new(cidx,:) = unique(mostLinearDep(cidx,:),'stable'); 
% end 
% mostLinearDep = mostLinearDep_new; 

dataproj(:,columnsToDelete) = []; 

    % Remake dataproj covariance matrix and ranking matrix 
    dCovMatrix(:,touchedColumns) = nan; 
    dCovMatrix(touchedColumns,:) = nan; 
    dCovMatrix(:,columnsToDelete) = []; 
    dCovMatrix(columnsToDelete,:) = []; 
    for bidx1 = 1:size(dataproj,2) 
        for bidx2 = (bidx1+1):size(dataproj,2)
            if(isnan(dCovMatrix(bidx1,bidx2)) ) 
                dCovMatrix(bidx1,bidx2) = det(dataproj(:,[bidx1 bidx2])'*dataproj(:,[bidx1 bidx2]));             
                dCovMatrix(bidx2,bidx1) = dCovMatrix(bidx1,bidx2); 
            end             
        end 
    end             
    bmax = max(max(abs(dCovMatrix))); 
    mostLinearDep = nan([1 1]*size(dCovMatrix,2)); % Each row is sorted, best candidate for xidx is at yidx=1 
    linearDepMeasure = nan([1 1]*size(dCovMatrix,2)); 
    for bidx = 1:size(dCovMatrix,2) 
        dCovMatrix(bidx,bidx) = 1.1*bmax; % Set diagonal to large value so we can find minimum of non-zero values 
        [linearDepMeasure(bidx,:),mld] = sort(dCovMatrix(bidx,:)); 
        mostLinearDep(bidx,:) = mld; 
    end 
        
nan(0); 
end

[v,d] = eig(qrproj*qrproj'); 
[d,sortfilter] = sort(diag(real(d)),'descend'); 
v = v(:,sortfilter); 



nan(0); 


