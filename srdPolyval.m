function [ polyTerms ] = srdPolyval( x, polyPowers )
% x: numDimensions -by- numElements 
% polyPowers: numTerms -by- numDimensions 

polyTerms = nan([size(x,2) size(polyPowers,1)]); 
for termidx = 1:size(polyPowers,1) 
    xcolumn_temp = ones([size(x,2) 1]); 
    for pcidx = 1:size(x,1) 
        if( 0 ~= polyPowers(termidx,pcidx) ) 
            xcolumn_temp = xcolumn_temp .* x(pcidx,:)'.^polyPowers(termidx,pcidx); 
        end 
    end     
    polyTerms(:,termidx) = xcolumn_temp; 
end


end

