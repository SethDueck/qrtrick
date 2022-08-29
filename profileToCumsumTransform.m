function [ transData ] = profileToCumsumTransform( data, numDatapoints )

cumsumData = cumsum(data,2); 

xq = linspace(min(min(cumsumData)),max(max(cumsumData)),numDatapoints+2); 
xq = xq(2:end-1); 

transData = nan([size(data,1),numDatapoints]); 
indexList = 1:size(cumsumData,2); 
interpolator = griddedInterpolant(indexList,indexList,'linear','none'); 
for didx = 1:size(data,1) 
    if(4343==didx)
        nan(0); 
    end 
    x = cumsumData(didx,:); 
%     xq = linspace(cumsumData(didx,1),cumsumData(didx,end),numDatapoints+2); 
    interpolator.GridVectors = {x};
%     transData(didx,:) = interp1(x,indexList,xq(2:end-1)); 
    transData(didx,:) = interpolator(xq); 
end
transData(isnan(transData)) = 0; 

end

