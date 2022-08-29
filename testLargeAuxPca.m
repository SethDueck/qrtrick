function [ result ] = testLargeAuxPca( x0,y0 )

nx = 100; 
ny = 100; 
zval = 1e0; 

data = [linspace(0,1,nx) linspace(0,0,ny); linspace(0,0,nx) linspace(0,1,ny)]';
data(:,1) = data(:,1) + x0 + 1e6; 
data(:,2) = data(:,2) + y0 + 1e6; 
[v,d] = eig(data'*data); 

data2 = nan([size(data,1) size(data,2)+1]); 
data2(:,1:2) = data; 
data2(:,3) = zval; 
[v2,d2] = eig(data2'*data2); 

result = [v2(1:2,end)*zval v2(end,1:2)'*zval]; 

end

