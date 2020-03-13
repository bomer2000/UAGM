% Uncertainty update and integration based on contact relationship

% if exist('pZlDeK')
    clear x m pZlD pZlDK pKlZ pZlDe pZlDeK pZlDeK_norm pZlDeK_cdf;
% end

Z_Pfieldmin= ;	%range of Z
Z_Pfieldmax= ;
dz_Pfield= ;    %resolution in Z
n_Pfield=ceil(Z_Pfieldmax-Z_Pfieldmin)/dz_Pfield; 

Interface_0_domain=zeros(ny_grid,nx_grid);
Interface_0_domain(find(~isnan(model_interface{1,1})==1))=1;
Interface_0_domain(1:50,:)=1;
%update of probability
for i=1:ny_grid
 for j=1:nx_grid
     for k=1:n_interface    
X_query=X0+(j-1)*d_grid;	
Y_query=Y0+(i-1)*d_grid;
if ~isnan(model_interface{k,1}(i,j))
zm=model_interface{k,1}(i,j);   %existing model
zkrig=InterfaceZ_xyKRIGU{k,1}(i,j);    %kriging prediction
zv=vInterfaceZ_xyKRIGU{k,1}(i,j);    %kriging prediction variance
if isnan(zkrig)
    zkrig=m_Interface(k,1);
    zv=v_Interface(k,1)^2;
end

for m=1:n_Pfield    
    ZSpace_Pfield(m,1)=Z_Pfieldmin+dz_Pfield*(m-1);   
    pZlD(m,1)=normpdf(ZSpace_Pfield(m,1),zkrig,sqrt(zv));   
    pZlDK(m,1)=normpdf(ZSpace_Pfield(m,1),zm,sqrt(zv));          
end

for m=1:n_Pfield    
    pKlZ(m,1)=exp(((zm-zkrig)*ZSpace_Pfield(m,1)+(zkrig^2-zm^2)/2)/zv);
end

% update BMEpdf with existing model
row=(j-1)*ny_grid+i;    
if isnan(InterfaceZ_xypdf{k,1}(row,1))  
    pZlDe=pZlD;
else    
    pZlDe=InterfaceZ_xypdf{k,1}(row,:)';
    pZlDe(isnan(pZlDe))=0;  
    pZlDe(find(pZlDe<0))=0; 
    pZlDe(find(pZlDe>1))=1; 
end
% show P(Z|X,E)
% plot(ZSpace_Pfield,pZlDe,'b:');
% hold on;
pZlDeK=pZlDe.*pKlZ;
pZlDeK(isnan(pZlDeK))=0;

c=1/trapz(ZSpace_Pfield,pZlDeK);
pZlDeK_norm=c*pZlDeK;   
pZlDeK_Pdf{i,j}(:,k)=pZlDeK_norm;
% show P(Z|X,E,M)
% plot(ZSpace_Pfield,pZlDeK_norm,'r--');
% title('P(Z|X)(r), P(Z|X,M)(g), P(Z|X,E)(b), P(Z|X,E,M)(r--)');

%transform pdf to cdf
for m=1:size(pZlDeK_norm,1)
    if m>1
        pZlDeK_cdf(m,1)=trapz(ZSpace_Pfield(1:m,1),pZlDeK_norm(1:m,1));
    end
end
% figure;
% plot(ZSpace_Pfield,pZlDeK_cdf,'b-');
% title('P(Z|X,E,M) Cdf');

pZlDeK_Pfieldcdf{i,j}(:,k)=pZlDeK_cdf;

Interface_ZPfield{k,1}(i,j)=ZSpace_Pfield(No_ZPfield,1);
else
    pZlDeK_Pdf{i,j}(:,k)=zeros(size(ZSpace_Pfield));
    pZlDeK_Pfieldcdf{i,j}(:,k)=zeros(size(ZSpace_Pfield));
end
     end
     pZlDeK_Pdf{i,j}(find(pZlDeK_Pdf{i,j}<0))=0;
     pZlDeK_Pfieldcdf{i,j}(find(pZlDeK_Pfieldcdf{i,j}<0))=0;
     pZlDeK_Pfieldcdf{i,j}(find(pZlDeK_Pfieldcdf{i,j}>1))=1;
     pZlDeK_Pfieldcdf{i,j}(isnan(pZlDeK_Pfieldcdf{i,j}))=0;

%Contact relationship of strata
Contact_list(1,1)=Interface_0_domain(i,j);	%starting from stratum ID: 0
for k=1:n_interface
    Contact_list(k+1,1)=~isnan(model_interface{k,1}(i,j));
end
k_add_1=1;    
ZPfield{i,j}(:,k_add_1)=ones(size(ZSpace_Pfield));

%Integration based on contact relationship
for k_add_1=2:n_interface+1    
    if Contact_list(k_add_1,1)  
        ZPfield{i,j}(:,k_add_1)=pZlDeK_Pfieldcdf{i,j}(:,k_add_1-1);
        if Contact_list(k_add_1-1,1)  
            ZPfield{i,j}(:,k_add_1-1)=ZPfield{i,j}(:,k_add_1-1).*(ones(size(ZSpace_Pfield))-ZPfield{i,j}(:,k_add_1));
            if k_add_1>2
                for k_2layer_before=1:k_add_1-2
                    temp_P=ZPfield{i,j}(:,k_add_1).*ZPfield{i,j}(:,k_2layer_before);
                end
                ZPfield{i,j}(:,k_add_1)=ZPfield{i,j}(:,k_add_1)-temp_P;
            end
        else    
            for k_1layer_before=1:k_add_1-1
                ZPfield{i,j}(:,k_1layer_before)=ZPfield{i,j}(:,k_1layer_before).*(ones(size(ZSpace_Pfield))-ZPfield{i,j}(:,k_add_1));
            end
        end
    else
        ZPfield{i,j}(:,k_add_1)=0;    
    end
end
    ZPfield{i,j}(find(ZPfield{i,j}<0))=0;
    ZPfield{i,j}(find(ZPfield{i,j}>1))=1;
 end
end
fprintf('Probability Integration Complete\n');
% show
% figure;
% plot(ZSpace_Pfield,ZPfield{i,j}(:,3),'r-');
% hold on;
% plot(ZSpace_Pfield,ZPfield{i,j}(:,4),'g-');
% hold on;
% plot(ZSpace_Pfield,ZPfield{i,j}(:,6),'b-');
% hold on;
% plot(ZSpace_Pfield,ZPfield{i,j}(:,8),'m-');
% hold on;
% plot(ZSpace_Pfield,ZPfield{i,j}(:,10),'c-');
