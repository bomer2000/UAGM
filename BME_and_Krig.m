%parameter setup before prediction(see details in BMElib)
X0= ;
Y0= ;
minc=[X0 Y0];            
dc=[ ];               
nc=[ ];            
[ck]=creategrid(minc,dc,nc);
nhmax=;
dmax=;
order=; 

%Error of drill data
error_drill= ;
%Error of section data
error_section= ;

%Transform determinate data to uncertain data(pdf data)
for i=1:n_interface   
    v_ZME{i,1}=error_section^2*ones(size(Section_point{i,1},1),1);   
    v_ZME{i,2}=error_guessdata(i,1)^2*ones(size(guess_data{i,1},1),1);  
    v_ZME{i,3}=error_drill^2*ones(size(Drill_point{i,1},1),1);  
	[softpdftype{i,1},nl{i,1},limi{i,1},probdens{i,1}]=probaGaussian(Section_point{i,1}(:,3),v_ZME{i,1});
    if size(guess_data{i,1},1)>0
        [softpdftype{i,2},nl{i,2},limi{i,2},probdens{i,2}]=probaGaussian(guess_data{i,1}(:,3),v_ZME{i,2}); 
	end
    [softpdftype{i,3},nl{i,3},limi{i,3},probdens{i,3}]=probaGaussian(Drill_point{i,1}(:,3),v_ZME{i,3}); 
end
%kriging (see details in BMElib)
for i=1:n_interface
    if size(guess_data{i,1},1)>0
    [InterfaceZ_xyKRIGU{i,1},vInterfaceZ_xyKRIGU{i,1}]=kriging(ck,[Interface_point{i,1}(:,1:2);guess_data{i,1}(:,1:2)],[Interface_point{i,1}(:,3);guess_data{i,1}(:,3)],modelInterfaceZ_xy_Cov{i},paramInterfaceZ_xy_Cov{i},nhmax,dmax,order);
    [ck1,ck2,InterfaceZ_xyKRIGU{i,1}]=col2mat(ck,InterfaceZ_xyKRIGU{i,1});
    [ck1,ck2,vInterfaceZ_xyKRIGU{i,1}]=col2mat(ck,vInterfaceZ_xyKRIGU{i,1});
    else
    [InterfaceZ_xyKRIGU{i,1},vInterfaceZ_xyKRIGU{i,1}]=kriging(ck,[Interface_point{i,1}(:,1:2)],[Interface_point{i,1}(:,3)],modelInterfaceZ_xy_Cov{i},paramInterfaceZ_xy_Cov{i},nhmax,dmax,order);
    [ck1,ck2,InterfaceZ_xyKRIGU{i,1}]=col2mat(ck,InterfaceZ_xyKRIGU{i,1});
    [ck1,ck2,vInterfaceZ_xyKRIGU{i,1}]=col2mat(ck,vInterfaceZ_xyKRIGU{i,1});
    end
end

nsmaxZ= ;	%maximum amount of softdata in neighbourhood
nhmaxZ= ;	%%maximum amount of harddata in neighbourhood
% BME prediction(BMEprobaPdf,see details in BMElib)
dmax_BME= ;
for i=1:n_interface 
    time=datestr(now,31);
    fprintf(' BME  Layer %d start: %s\n',i,time);
if size(guess_data{i,1},1)>0
    for j=1:size(ck,1)
        [InterfaceZ_xyBMEprobaPdf{i,1}(j,:),InterfaceZ_xypdf{i,1}(j,:),InterfaceZ_xyinfo{i,1}(j,:)]=BMEprobaPdf(ZSpace_Pfield,ck(j,1:2),Drill_point{i,1}(:,1:2),[Section_point{i,1}(:,1:2);guess_data{i,1}(:,1:2)],Drill_point{i,1}(:,3),[softpdftype{i,1}],[nl{i,1};nl{i,2}],[limi{i,1};limi{i,2}],[probdens{i,1};probdens{i,2}],modelInterfaceZ_xy_Cov{i},paramInterfaceZ_xy_Cov{i},nhmaxZ,nsmaxZ,dmax_BME,order);
        continue
    end
else
    for j=1:size(ck,1)
        [InterfaceZ_xyBMEprobaPdf{i,1}(j,:),InterfaceZ_xypdf{i,1}(j,:),InterfaceZ_xyinfo{i,1}(j,:)]=BMEprobaPdf(ZSpace_Pfield,ck(j,1:2),Drill_point{i,1}(:,1:2),Section_point{i,1}(:,1:2),Drill_point{i,1}(:,3),[softpdftype{i,1}],[nl{i,1}],[limi{i,1}],[probdens{i,1}],modelInterfaceZ_xy_Cov{i},paramInterfaceZ_xy_Cov{i},nhmaxZ,nsmaxZ,dmax_BME,order);
        continue
    end
end
time=datestr(now,31);
fprintf('BME  Layer %d complete: %s\n',i,time); 
end
