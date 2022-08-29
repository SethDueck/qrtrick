
for cidx = 1:numel(cindices_fits) 
    if(0==numel(cindics_fits{cidx})) 
        continue; 
    end 
    
sd = ssdata(cindices_fits{cidx}(f),:); 
krigPoints = randperm(size(sd,1),50); 
% krigPoints = nan([numel(cindices_final) 1]); 
% for nidx = 1:numel(cindices_final) 
%     krigPoints(nidx) = find(filter_worst==cindices_final{nidx}(end)); 
% end 

khr = hr(cindices_fits{cidx}(c5f),40); 
kinterpvalues = khr(krigPoints); 

temp = squareform(pdist(sd(krigPoints,:))); 

% A*lambda = b
kScaleFac = -5/max(temp(:));
kA = nan([1 1]*(numel(krigPoints)+1));
kA(1:numel(krigPoints), 1:numel(krigPoints)) = exp(kScaleFac*temp); 
kA(end,:) = 1; 
kA(:,end) = 1; 
kA(end,end) = 0; 
invkA = inv(kA); 


kb = nan([numel(krigPoints)+1 1]); 
kb(end) = 1; 
kresult = nan(size(khr)); 
kLambdas_all = nan([numel(krigPoints) size(sd,1)]); 
for pidx = 1:size(sd,1) 
    temp = pdist2(sd(krigPoints,:),sd(pidx,:)); 
    kb(1:end-1) = exp(kScaleFac*temp); 
    kLambdas = invkA*kb; 
    kresult(pidx) = sum(kLambdas(1:end-1).*kinterpvalues); 
    kLambdas_all(:,pidx) = kLambdas(1:end-1); 
end 

end 



