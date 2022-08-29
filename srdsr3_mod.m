function [ xk, wk ] = srdsr3_mod( A, C, b, pointWeights, eta, lambda, w0, convergenceThreshold, maxNumIterations )

Apw = A; 
if( any(1~=pointWeights) ) 
    % Only do this step if necessary; it can actually take the bulk of
    % runtime if there are many short calls to this function. 
    for aidx = 1:size(A,1) 
        Apw(aidx,:) = Apw(aidx,:)*pointWeights(aidx); 
    end 
end 

% Atb = A'*(pointWeights.*b); 
Atb = Apw'*b; 
Hinv = inv(A'*Apw + (1/eta)*(C'*C)); 

F_kappa = [(1/eta)*Apw*(Hinv*C'); sqrt(1/eta)*(eye(size(C,1))-(1/eta)*C*Hinv*C')]; 
FtF = F_kappa'*F_kappa; 
% G_kappa = [diag(pointWeights)*(eye(size(A,1)) - A*Hinv*A'); sqrt(1/eta)*C*Hinv*A'] * diag(pointWeights); 
% g_kappa = G_kappa*b; 
% Equivalent to lines above but avoids an nPoints*nPoints memory
% allocation: 
% g_kappa =  [ pointWeights.* (pointWeights.*b - A*(Hinv*(A'*(pointWeights.*b)))) ; ... 
%     sqrt(1/eta)*C*(Hinv*(A' * (pointWeights.*b)))]; 
g_kappa =  [ pointWeights.* (pointWeights.*b - A*(Hinv*(Atb))) ; ... 
    sqrt(1/eta)*C*(Hinv*Atb)]; 
Ftg_kappa = F_kappa'*g_kappa; 


k = 0; 
wk = w0; 
notConverged = true; 
allwk = nan([numel(w0) maxNumIterations]); 
while notConverged
    k = k+1; 

    wk_old = wk; 
%     deltaWk = - eta*F_kappa'*(F_kappa*wk_old - g_kappa); 
%     temp = (F_kappa*wk_old - g_kappa); 
%     deltaWk = - eta*F_kappa'*temp; 
    deltaWk = -eta*(FtF*wk_old) + eta*Ftg_kappa; 
    
    proxArg = wk_old + deltaWk; 
    wk = sign(proxArg).*max([abs(proxArg)-lambda,zeros(size(proxArg))],[],2); % Soft thresholding  
    allwk(:,k) = wk; 
    
    sumNetDeltaWk = sum(abs(wk-wk_old)); 
    notConverged = (convergenceThreshold < sumNetDeltaWk) && (k<maxNumIterations); 
end 

xk = (Hinv) * (Atb + (1/eta)*C'*wk); 


nan(0); 

end

