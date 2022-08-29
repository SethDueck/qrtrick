function [ u,svec,v ] = srd_svd( data, cutoffFactor )
% Returns u*diag(svec)*v' = data; 

if( size(data,1) < size(data,2) )
    data = data'; 
end 

[v,d] = eig(data'*data); 
[d,sortfilter] = sort(real(diag(d)), 'descend'); 
sortfilter(d<=sum(d)*cutoffFactor) = []; 
d(d<=sum(d)*cutoffFactor) = []; 
v = v(:,sortfilter); 

svec = sqrt(d); 
temp = v*diag(1./svec); 
u = data*(temp); 

end

