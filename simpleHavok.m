

% Make dummy dataset 
t = linspace(0,1,1e4); 
omega1 = 80; 
omega2 = 100; 
% y_real = 0.5*t + 0.8*t.*t + 2; 
% % y_osc = 0.05*cos(omega1*t); 
% y_osc = 0.05*cos(omega1*t + omega2*t.*t); 
% y_noise = (max(y_real)-min(y_real))*randn(size(y_real)) * 1e-3; % Add noise 
y_real = atan(21*t-4) + atan(48*t-34); 
y_osc = cos(omega1*t + omega2*t.*t) * 1e-1; 
y_noise = (max(y_real)-min(y_real))*randn(size(y_real)) * 1e-3; % Add noise 
y = y_real + y_osc + y_noise; 

% Construct hankel matrix, there's probably a faster way to do this  
hankelWindowSize = 1e3; 
h = nan([hankelWindowSize (numel(y)-hankelWindowSize+1)]); 
for hidx = 1:(numel(y)-hankelWindowSize+1) 
    h(:,hidx) = y(hidx-1+(1:hankelWindowSize));
%     h(:,hidx) = h(:,hidx) - mean(h(:,hidx)); 
end 

% Do taken's embedding, this is equivalent to what he does in the video but
% is faster than calling svd/svds: 
[u,s] = eigs(h*h',20); 
s = sqrt(s); 
v = h'*(u*diag(1./diag(s))); 
numImportantDims = sum( (diag(s.^2)/sum(diag(s).^2)) > 1e-6 ); 

% Take dv/dt, average v for each timestep (I get to assume I have constant
% time steps, make the obvious adjustment if that's not true for you) 
dvdt = diff(v(:,1:numImportantDims)); 
vmean = ( v(1:end-1,1:numImportantDims) + v(2:end,1:numImportantDims) ) / 2; 

% Do sindy for each coordinate. Here I break a little with the notation he
% uses in the video because he uses like 'A' in a few different places 
regressionModel = nan([1 1]*size(vmean,2)); 
pointWeights = ones([size(dvdt,1) 1]); 
for vidx = 1:size(regressionModel,1) 
    A = vmean; 
    b = dvdt(:,vidx); 
    C = eye(size(regressionModel,1)); 
    
    eta = 1e0; 
    lambda = 1e-3; 
    w0 = ones([size(A,2) 1]); 
    convergenceThreshold = 1e-5; 
    maxNumIterations = 1e4; 
    [ xk, wk ] = srdsr3_mod( A, C, b, pointWeights, eta, lambda, w0, convergenceThreshold, maxNumIterations ); 
    
    regressionModel(vidx,:) = wk; 
    
    nan(0); 
    
end 
 

% Look for rotations in the regression model 
[v_model,d_model] = eig(regressionModel);
d_model = diag(d_model); 
[~,sortfilter] = sort(abs(imag(d_model))); % Sort from slow-oscillating to fast-oscillating 
d_model = d_model(sortfilter); 
v_model = v_model(:,sortfilter); 
d_conjPartnerIdx = nan(size(d_model)); 
d_conjPairK = nan(size(d_model)); 
dcpkidx = 1; 
for didx = 1:numel(d_model) 
    if(abs(imag(d_model(didx)))>1e-10) 
        % This eigenvector has a rotation pair 
        [~,d_conjPartnerIdx(didx)] = min(abs(imag(d_model+d_model(didx)))); 
        if( isnan(d_conjPairK(d_conjPartnerIdx(didx))) ) 
            d_conjPairK(didx) = dcpkidx; 
            d_conjPairK(d_conjPartnerIdx(didx)) = dcpkidx; 
            dcpkidx = dcpkidx+1; 
        else 
            if(d_conjPairK(didx)~=d_conjPairK(d_conjPartnerIdx(didx))) 
                warning('I assumed each rotation subspace had dimension 2, looks like I was wrong -- code needs simple mods to handle more generalized eigenvectors.'); 
            end 
        end
    else 
        % This eigenvector is an axis of rotation, call all of these parts 
        % "group 0" 
        d_conjPartnerIdx(didx) = 0; 
        d_conjPairK(didx) = 0; 
    end 
end 
% Recast phase-space coordinates onto the eigenvectors of the koopman
% operator 
p_model = v_model'*v(:,1:numImportantDims)'; 
% Oscillating parts will store potential energy -- treat ode elements with
% low potential energy as part of the "smooth" signal (k=0) 
d_imaginaryEnergy = sum((real(diag(d_model)*p_model).*imag(diag(d_model)*p_model)).^2,2);    
maxEnergyFrac = 1e-4; 
filter_treatAsReal = (d_imaginaryEnergy/sum(d_imaginaryEnergy)) < maxEnergyFrac; 
[~,~,uniqueFrequenciesK] = unique(abs(imag(d_model))); % Remember how d_model was sorted near the top of this code block 
numLowFreqs0 = max(uniqueFrequenciesK(filter_treatAsReal)); 
d_conjPairK( numLowFreqs0 >= uniqueFrequenciesK ) = 0; 
d_conjPairK( numLowFreqs0 <  uniqueFrequenciesK ) = uniqueFrequenciesK( numLowFreqs0<uniqueFrequenciesK) - numLowFreqs0; 



% We grouped together parts that oscillate in the same subspace; convert
% each of these parts back to signal-space seperately. There's definitely a
% faster way to do this but I don't feel like working on it. 
y_ret_parts = nan([max(d_conjPairK)+1 ,numel(y)]); 
for kidx = 0:max(d_conjPairK) 
    p_model_k = p_model; 
    p_model_k( kidx~=d_conjPairK,:) = 0; 
    vThisK = v_model' \ p_model_k; 
    
    temp = flipud((u(:,1:numImportantDims)*s(1:numImportantDims,1:numImportantDims))*vThisK); 
    for yidx = 1:size(y_ret_parts,2) 
        y_ret_parts(kidx+1,yidx) = mean(diag(temp,-size(h,1)+yidx)); 
    end
end 

% The result is guaranteed to be real, but some imaginary stuff survives
% due to machine error 
y_ret_parts = real(y_ret_parts); 

% This is the same as how I did y_ret_real, just inverting from taken's 
% space to signal-space without worrying about eigenvectors of the koopman
% operator 
temp = flipud(u*s*v'); 
y_ret = nan(size(y)); 
for yidx = 1:numel(y_ret) 
    y_ret(yidx) = mean(diag(temp,-size(h,1)+yidx)); 
end 


%
% Plot results 
%
figure(1911) 
clf 
hold all 
plot(t,y,'k','linewidth',2.5) 
plot(t,y_ret_parts(1,:),'Color',[1 0.2 0.9]*0.9,'LineWidth',2.5) 
for kidx = 2:size(y_ret_parts,1) 
    plot(t,y_ret_parts(1,:)-y_ret_parts(kidx,:),'linewidth',1.5)
end 
plot(t,y,'k','linewidth',2.5) 
plot(t,y_ret_parts(1,:),'Color',[1 0.2 0.9]*0.9,'LineWidth',2.5) % Plot this twice for visibility, annoying 
title('Signal and retrieved signal') 
legend('Signal','Non-oscillating retrieved part',['All other lines are non-oscillating' char(10) 'part MINUS one oscillating part' char(10) 'of the retrieved signal'])


figure(1912) 
clf 
hold all 
plot(y_ret_parts') 
title('Each component of the retrieved signal') 


