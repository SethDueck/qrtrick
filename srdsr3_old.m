function [ xk, wk ] = srdsr3( A, C, b, eta, lambda, w0, convergenceThreshold, maxNumIterations )

notConverged = true;  

AtA = A'*A; 
CtC = C'*C; 
Atb = A'*b; 
invTerm = inv(AtA+CtC/eta); 

k = 0; 
wk = w0; 
while notConverged
    k = k+1; 
    xk = (invTerm)*(Atb+C'*wk); 
    
    % Prox operator 
    Cx = C*xk; 
    wk_old = wk; 
    wk = sign(Cx).*max([abs(Cx)-lambda,zeros(size(Cx))],[],2); % Soft thresholding
    notConverged = (convergenceThreshold < sum(abs(wk-wk_old))) && (k<maxNumIterations); 
end 

xk = (invTerm)*(Atb+C'*wk); 

end

