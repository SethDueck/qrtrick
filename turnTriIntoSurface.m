function [ s_unique ] = turnTriIntoSurface( t )
% t is output to conv hull, probably nHyperTriangles-by-4. 
% s is the list of all 3-vertex triangles in t 

% Sort indices of vertices in t; make s(1,X)<s(2,X) below 
t = sort(t,2); 

% Index all edges (pairs of vertices) in a hypertriangle 
smallLoopSize = nchoosek(size(t,2),3); 
smallLoopIndices = false([smallLoopSize size(t,2)]); 
sidx = 1; 
for vidx1 = 1:size(t,2) 
    for vidx2 = (vidx1+1):size(t,2) 
        for vidx3 = (vidx2+1):size(t,2) 
            smallLoopIndices(sidx,vidx1) = true; 
            smallLoopIndices(sidx,vidx2) = true; 
            smallLoopIndices(sidx,vidx3) = true; 
            sidx = sidx + 1;     
        end 
    end 
end 

% Move edges from t to s 
s = nan([3, size(t,1)*smallLoopSize]); 
sidx = 1; 
for tidx1 = 1:size(t,1) 
    for tidx2 = 1:smallLoopSize 
        s(:,sidx) = t(tidx1,smallLoopIndices(tidx2,:)); 
        sidx = sidx + 1; 
    end 
    
end 
nan(0); 

% Remove duplicate edges 
s_unique = unique(s','rows'); 


end

