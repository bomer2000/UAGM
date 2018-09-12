%filter the posterior probability field with DTM
ZPfield_DTM=ZPfield;
for i=1:ny_grid
 for j=1:nx_grid
     ZPfield_DTM{i,j}(find(ZSpace_Pfield>DTM(i,j)),:)=-1;
 end
end
fprintf('ZPfield DTM filtration Complete\n');
