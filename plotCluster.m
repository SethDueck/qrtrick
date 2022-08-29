

cidx = 1; 

for cidx = 1:218 
    
    thisCluster = ssdata(cidx==clabels,:); 
    tcn = thisCluster - repmat(mean(thisCluster,1),[size(thisCluster,1) 1]); 
    [vcn,dcn] = eig(tcn'*tcn); 
    [dcn,sortfilter] = sort(real(diag(dcn)),'descend'); 
    vcn = vcn(:,sortfilter); 
    pcn = vcn'*tcn'; 
    
    figure(10) 
    clf 
    
    subplot(1,2,1)
    surf(thisCluster) 
    title([num2str(cidx) ': ' num2str(size(thisCluster,1))]) 

    subplot(1,2,2) 
    plot3(pcn(3,:),pcn(2,:),pcn(1,:),'.'); 
    axis equal 
    
    pause() 
end 


