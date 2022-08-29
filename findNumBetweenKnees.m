function [ count ] = findNumBetweenKnees( data )

count = sum( (data<findKnee(data)) & (data>-findKnee(-data)) ); 

end

