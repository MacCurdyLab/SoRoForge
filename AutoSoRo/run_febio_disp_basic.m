function [V_DEF,timeVec] = run_febio_disp_basic(F,V,elementMaterialID,appliedPressure,FL)

fixed_set = find(V(:,3)<0.05*max(V(:,3)));
tip_set = find(V(:,3)>0.95*max(V(:,3)));
w_t = 2.2;
%appliedPressure=1e-2; 

%% Control parameters

% Path names
defaultFolder = fileparts(mfilename('fullpath'));
savePath=fullfile(defaultFolder,'data','temp');
% mkdir('data');mkdir('data/temp');
mkdir([defaultFolder '/data']);
mkdir([defaultFolder '/data/temp']);

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement

%Material parameter set
shoreA = 80;
E = 10^((shoreA)*0.0235-0.6403);
E = 3;
v = 0.49;
G = E/(2*(1+v));
c1=0.03;        %Shear-modulus-like parameter
c1=G;
m1=2;           %Material parameter setting degree of non-linearity
k_factor=1e2;   %Bulk modulus factor 
k=c1*k_factor;  %Bulk modulus

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=10; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=12;%12; %Optimum number of iterations
max_retries=8;%10; %Maximum number of retires
dtmin=(1/numTimeSteps)/4; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

F_pressure=FL;


%% Split the mesh into 2 groups
E = F;
if range(elementMaterialID) == 0
    elementMaterialID(1) = elementMaterialID(1)+1;
end
[elementMaterialID,I] = sort(elementMaterialID,'ascend');
E = E(I,:);

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_ups=max_ups;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

%Material section
materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c1;
febio_spec.Material.material{2}.m1=m1;
febio_spec.Material.material{2}.c2=c1;
febio_spec.Material.material{2}.m2=-m1;
febio_spec.Material.material{2}.k=k;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='SoftPart';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tri3'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:nnz(elementMaterialID==1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E(elementMaterialID==1,:);

partName2='StiffPart';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tri3'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=(1+nnz(elementMaterialID==1):1:nnz(elementMaterialID==1)+nnz(elementMaterialID==2))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E(elementMaterialID==2,:);

meshDataName1='ThicknessData1';           
febio_spec.MeshData.ElementData{1}.ATTR.name=meshDataName1;
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:nnz(elementMaterialID == 1))';
febio_spec.MeshData.ElementData{1}.elem.VAL=w_t*ones(nnz(elementMaterialID == 1),3);

meshDataName2='ThicknessData2';           
febio_spec.MeshData.ElementData{2}.ATTR.name=meshDataName2;
febio_spec.MeshData.ElementData{2}.ATTR.elem_set=partName2;
febio_spec.MeshData.ElementData{2}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{2}.elem.ATTR.lid=(1:nnz(elementMaterialID == 2))';
febio_spec.MeshData.ElementData{2}.elem.VAL=w_t*ones(nnz(elementMaterialID == 2),3);

% -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL=F_pressure;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=fixed_set(:);

nodeSetName2='tipDispList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=tip_set(:);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';

%Loads section
% -> Surface load    
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=-appliedPressure;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=0;
febio_spec.Loads.surface_load{1}.linear=0;
febio_spec.Loads.surface_load{1}.shell_bottom=0;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|
% febView(febio_spec)

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='external';%'internal';
febioAnalysis.t_check=5; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=20; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 
    
try 
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    
    timeVec=dataStruct.time; %Time
        
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    
catch
    V_DEF = [];
    timeVec = [];
end
