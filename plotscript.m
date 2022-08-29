


nidx = 6; 

figure(8) 
clf 
hold all 
% scatter3(Pall{nidx}(3,:),Pall{nidx}(2,:),Pall{nidx}(1,:),2,pscores_exp(nidx,:))
scatter3(Pall{nidx}(3,:),Pall{nidx}(2,:),Pall{nidx}(1,:),2,max(pscores_exp,[],1))
axis vis3d


