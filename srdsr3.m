function [ xk, wk ] = srdsr3( A, C, b, eta, lambda, w0, convergenceThreshold, maxNumIterations )


Atb = A'*b; 
Hinv = inv(A'*A + (1/eta)*(C'*C)); 

F_kappa = [(1/eta)*A*Hinv*C'; sqrt(1/eta)*(eye(size(C,1))-(1/eta)*C*Hinv*C')]; 
G_kappa = [eye(size(A,1)) - A*Hinv*A'; sqrt(1/eta)*C*Hinv*A']; 
g_kappa = G_kappa*b; 

k = 0; 
wk = w0; 
notConverged = true; 
allwk = nan([numel(w0) maxNumIterations]); 
while notConverged
    k = k+1; 

    wk_old = wk; 
    deltaWk = - eta*F_kappa'*(F_kappa*wk_old - g_kappa); 

    proxArg = wk_old + deltaWk; 
    wk = sign(proxArg).*max([abs(proxArg)-lambda,zeros(size(proxArg))],[],2); % Soft thresholding  
    allwk(:,k) = wk; 
    
    sumNetDeltaWk = sum(abs(wk-wk_old)); 
    notConverged = (convergenceThreshold < sumNetDeltaWk) && (k<maxNumIterations); 
end 

xk = (Hinv) * (Atb + (1/eta)*C'*wk); 


nan(0); 

end

