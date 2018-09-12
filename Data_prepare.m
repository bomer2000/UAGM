%Input data and calculate covariance function
%Details about covariance calculation functions could be found in BMElib(http://www.unc.edu/depts/case/BMELIB/).
%Input Interface point data(vector of X,Y,Z,stratum ID,source(1:drill£¬2:section))
Interface_data= ;
%Input drill data(vector of X,Y,Z,stratum ID,source(1:drill£¬2:section),drill ID)
Drill_data= ;
%Input section data(vector of X,Y,Z,stratum ID,source(1:drill£¬2:section),section ID)
Section_data= ;
error_section= ;    %error of section data
%Read DTM data(X,Y,Z)
DTM_data= ;
%Read  interface model
model_interface{stratum ID,1}= ;%raster of elevation value of the geological interface 

%Auxiliary data form experts' guess
guess_data{stratum ID,1}=[ , , ];%(X,Y,Z)

n_data=size(Interface_data,1);  %amount of Interface points
n_drill=size(Drill_data,1);   %amount of Drill_data
n_section=size(Section_data,1);   %amount of Section_data
n_interface= ; %amount of Interfaces
n_layer=n_interface+1;    %amount of strata

%Error of guess_data
for i=1:n_interface	
    error_guessdata(i,1)= ;
end

%Study area setup
Xmin= ;
Xmax= ;
Ymin= ;
Ymax= ;
Zmin= ;
Zmax= ;

%Data classification
for i=1:n_interface
    Interface_point{i,1}=Interface_data(find(Interface_data(:,4)==i),1:3);  
    Drill_point{i,1}=Drill_data(find(Drill_data(:,4)==i),1:3);
    Section_point{i,1}=Section_data(find(Section_data(:,4)==i),1:3);
end
%Statistical information
for i=1:n_interface
    if size(guess_data{i,1},1)>0
        m_Interface(i,1)=mean([Interface_point{i,1}(:,3);guess_data{i,1}(:,3)]);
        v_Interface(i,1)=std([Interface_point{i,1}(:,3);guess_data{i,1}(:,3)]);
        s_Interface(i,1)=skewness([Interface_point{i,1}(:,3);guess_data{i,1}(:,3)]);
        k_Interface(i,1)=kurtosis([Interface_point{i,1}(:,3);guess_data{i,1}(:,3)]);
    else
        m_Interface(i,1)=mean([Interface_point{i,1}(:,3)]);
        v_Interface(i,1)=std([Interface_point{i,1}(:,3)]);
        s_Interface(i,1)=skewness([Interface_point{i,1}(:,3)]);
        k_Interface(i,1)=kurtosis([Interface_point{i,1}(:,3)]);
    end
end
for i=1:n_interface
        m_model(i,1)=mean([model_interface{i,1}(find(~isnan(model_interface{i,1})))]);
        v_model(i,1)=std([model_interface{i,1}(find(~isnan(model_interface{i,1})))]);
        s_model(i,1)=skewness([model_interface{i,1}(find(~isnan(model_interface{i,1})))]);
        k_model(i,1)=kurtosis([model_interface{i,1}(find(~isnan(model_interface{i,1})))]);
end
%the farthest distance of data point
for i=1:n_interface
    dxymax(i,1)=((max(Interface_point{i,1}(:,1))-min(Interface_point{i,1}(:,1)))^2+(max(Interface_point{i,1}(:,2))-min(Interface_point{i,1}(:,2)))^2)^(1/2);
    dzmax(i,1)=max(Interface_point{i,1}(:,3))-min(Interface_point{i,1}(:,3));
end

%Mian direction and secondary direction anisotropy
c_Mian=( : : )'; 	%vector giving the limits of the distance classes that are used for estimating the covariance function (see details in BMElib).
c_Second=( : : )';
options_Mian=[0, , ];  % Angles of anisotropy (see details in BMElib)
options_Second=[0, , ];

%Covariance calculation function: covario (see details in BMElib)
for i=1:n_interface
    if size(guess_data{i,1},1)>0
    [dInterfaceZ_Mian_Cov{i,1},vInterfaceZ_Mian_Cov{i,1},oInterfaceZ_Mian_Cov{i,1}]=covario([Interface_point{i,1}(:,1:2);guess_data{i,1}(:,1:2)],[Interface_point{i,1}(:,3);guess_data{i,1}(:,3)],c_Mian,'loop',options_Mian);
    [dInterfaceZ_Second_Cov{i,1},vInterfaceZ_Second_Cov{i,1},oInterfaceZ_Second_Cov{i,1}]=covario([Interface_point{i,1}(:,1:2);guess_data{i,1}(:,1:2)],[Interface_point{i,1}(:,3);guess_data{i,1}(:,3)],c_Second,'loop',options_Second);
    else
    [dInterfaceZ_Mian_Cov{i,1},vInterfaceZ_Mian_Cov{i,1},oInterfaceZ_Mian_Cov{i,1}]=covario(Interface_point{i,1}(:,1:2),Interface_point{i,1}(:,3),c_Mian,'loop',options_Mian);
    [dInterfaceZ_Second_Cov{i,1},vInterfaceZ_Second_Cov{i,1},oInterfaceZ_Second_Cov{i,1}]=covario(Interface_point{i,1}(:,1:2),Interface_point{i,1}(:,3),c_Second,'loop',options_Second);
    end
end

%Covariance function setup (see details in BMElib)
modelInterfaceZ_xy_Cov{interface ID}={'nuggetC','gaussianC'};
paramInterfaceZ_xy_Cov{interface ID}={[ ],[ ]};
modelplot(dInterfaceZ_xy_Cov{interface ID},modelInterfaceZ_xy_Cov{interface ID},paramInterfaceZ_xy_Cov{interface ID});

