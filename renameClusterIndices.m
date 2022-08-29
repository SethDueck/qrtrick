function [ clabels ] = renameClusterIndices(clabels, subclusterIndices, cidx) 

% subclusterIndices==1 remains the same. 
% All others get moved to the end of the clabel order 
subclusterIndices(subclusterIndices>1) = max(clabels) + subclusterIndices(subclusterIndices>1) - 1; 
subclusterIndices(subclusterIndices==1) = cidx; 

% Assign new labels 
clabels(clabels==cidx) = subclusterIndices; 

end

