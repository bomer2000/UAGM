%calculate the entropy field with ZPfield_DTM
for i=1:ny_grid
 for j=1:nx_grid
     Entropy_field_DTM{i,j}=zeros(n_Pfield,1);
     for m=1:n_Pfield   
         if ZPfield_DTM{i,j}(m,1)==-1
             Entropy_field_DTM{i,j}(m,1)=-1;
         else
             temp_entropy=ZPfield_DTM{i,j}(m,:).*log2(ZPfield_DTM{i,j}(m,:));
             temp_entropy(isnan(temp_entropy))=0;   
             Entropy_field_DTM{i,j}(m,1)=-sum(temp_entropy);
         end
     end
 end
end
fprintf('Entropy_field_DTM Complete\n');
