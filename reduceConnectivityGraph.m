function [ D ] = reduceConnectivityGraph( C )
% C is nClusters-by-nClusters. C(aidx,bidx) is the minimum distance between
% any point in cluster aidx to any point in cluster bidx. 
% Assume movement is free within clusters; output is a connectivity graph
% with the shortest distances between clusters (i.e. elements are the sum
% of "inter-cluster jumps" between clusters). 

D = C; 

isConverged = false; 


while ~isConverged 
    isConverged = true; 
    
    for didx1 = 1:size(D,1) 
        for didx2 = (didx1+1):size(D,2) 
            % Minimum distance with jump. If you first jump from didx1 to 
            % didx2, or vice versa, this is the minimum distance to all 
            % other clusters. 
            mdwj = D(didx1,didx2) + min(D(:,didx1),D(:,didx2)); 
            
            % Replace "direct" distance with mdwj if using the jump is
            % shorter
            for didx3 = 1:size(D,1) 
                if(     mdwj(didx3) < D(didx3,didx1) ) 
                    D(didx3,didx1) = mdwj(didx3); 
                    D(didx1,didx3) = mdwj(didx3); 
                    isConverged    = false; 
                elseif( mdwj(didx3) < D(didx3,didx2) ) 
                    D(didx3,didx2) = mdwj(didx3); 
                    D(didx2,didx3) = mdwj(didx3); 
                    isConverged    = false; 
                end 
            end
        end 
    end 
    
end 


end

