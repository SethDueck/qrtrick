

numQRVecs = 20; 
qrdata = ssdata; 

qrvecs   = nan(size(qrdata,1),numQRVecs); 
qrproj = nan(numQRVecs,size(qrdata,2)); 
qrIndices = nan([numQRVecs 1]); 

for nidx = 1:numQRVecs 
    [~,maxidx] = max(sum(abs(qrdata),1)); 
    qrIndices(nidx) = maxidx; 
    
    qrunit = qrdata(:,maxidx); 
    qrunit = qrunit / sqrt(sum(qrunit.^2)); 
%     qrscores = qrdata*qrunit'; 
    qrscores = qrunit'*qrdata; 
    
    qrvecs(:,nidx) = qrunit; 
    qrproj(nidx,:) = qrscores; 
    qrdata = qrdata - qrunit*qrscores; 
    
    nan(0); 
    
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
dataproj_orig = dataproj; % Copy original state of dataproj 
% TODO should replace lines below with #qrtrick_makeDataCovMatrix 
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
target_numDims = 3; 
linidx = 2; 
while size(dataproj,2)>target_numDims 
    if( size(dataproj,2)>(linidx+1) ) 
        linidx = linidx+1; 
    else 
        linidx = 2; 
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
end 

if( any( ismember(columnsToDelete,touchedColumns) ) ) 
    warning('Replacements and to-delete overlaps') 
end 
if( numel(columnsToDelete) ~= numel(unique(columnsToDelete)) ) 
    warning('Columns deleted more than once') 
end 

dataproj(:,columnsToDelete) = []; 

    % Remake dataproj covariance matrix and ranking matrix 
    dCovMatrix(:,touchedColumns) = nan; 
    dCovMatrix(touchedColumns,:) = nan; 
    dCovMatrix(:,columnsToDelete) = []; 
    dCovMatrix(columnsToDelete,:) = []; 
    
    [ temp1, temp2] = qrtrick_makeDataCovMatrix( dataproj, dCovMatrix ); 
    
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

% Find sensor points that best capture behaviour of original dataset 
% [v,d] = eig(dorig'*dorig); 
% [d,sortfilter] = sort(real(diag(d)),'descend'); 
% v = v(:,sortfilter); 
sensorPlacementData = dataproj_orig; 
numPivotPoints = 50; 
if(numPivotPoints>size(sensorPlacementData,2)) 
    % Cap number of sensors at number of dimensions in placement data; 
    % don't fit to noise. 
    numPivotPoints = size(sensorPlacementData,2); 
end 
pivotPointList = nan([numPivotPoints 1]); 
for pidx = 1:numPivotPoints 
    [~,midx] = max(sum(abs(sensorPlacementData),2)); 
    pivotPointList(pidx) = midx; 
    pivotUnit = sensorPlacementData(midx,:); 
    pivotUnit = pivotUnit / sqrt(sum(pivotUnit.*pivotUnit)); 
    pivotScore = sensorPlacementData*pivotUnit'; 
    sensorPlacementData = sensorPlacementData - pivotScore*pivotUnit; 
end 

% Try to make sharp/flattened auxiliary vectors 
for didx = 1:size(dataproj_orig,2) 
    for tempidx = 1:5
        exp(-dataproj_orig); 
    end
end 

% Make basis for decomposing radiances 
[vret,dret] = eig(dataproj_orig'*dataproj_orig); 
[dret,sortfilter] = sort(real(diag(dret)),'descend'); 
vret = vret(:,sortfilter); 
vp = dataproj_orig*vret; 
auxv = [ ...
      linspace( 1, 1, size(vp,1))', ... 
      linspace( 0, 1, size(vp,1))', ... 
      linspace( 0, 1, size(vp,1))'.^2 ... 
%       linspace(-1, 1, size(vp,1))'.^2 ... 
%       linspace(-1, 0, size(vp,1))'.^2 ... 
      ]; 
%   auxv = [];
vp = [auxv, vp]; 
vp_small = vp(pivotPointList,:); 
for didx = 1:size(vp,2)  
    vp(:,didx) = vp(:,didx) / sqrt(sum(vp(:,didx).^2)); 
    vp_small(:,didx) = vp_small(:,didx) / sqrt(sum(vp_small(:,didx).^2)); 
    for didx0 = 1:didx-1 % re-orthogonalize data after adding in aux dims 
        vp(:,didx) = vp(:,didx) - vp(:,didx0)*(vp(:,didx)'*vp(:,didx0));
        vp_small(:,didx) = vp_small(:,didx) - vp_small(:,didx0)*(vp_small(:,didx)'*vp_small(:,didx0)); 
    end 
    vp(:,didx) = vp(:,didx) / sqrt(sum(vp(:,didx).^2)); 
    vp_small(:,didx) = vp_small(:,didx) / sqrt(sum(vp_small(:,didx).^2));
end 

sr3basis = [dataproj,ones([size(dataproj,1) 1])]; 

% fakerad = dataproj_orig(:,4); 
fakerad = sum(vp,2); 
fakerad = vp(:,4); 
fakerad = fs(:,70)-ss(:,70); 
xk = qrtrick_invert(vp_small, fakerad(pivotPointList,:)); 
xk = qrtrick_invert(sr3basis(pivotPointList,:),fakerad(pivotPointList,:)); 

figure(8) 
clf 
hold all 
% scatter3(dataproj(:,3),dataproj(:,2),dataproj(:,1),2,fakerad) 
plot(fakerad)
plot((fakerad'*vp)*vp') 
% plot((fakerad(pivotPointList,:)'*vp_small)*vp')
plot(sr3basis*xk); 

nan(0); 


