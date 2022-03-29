function Params = setParams(varargin)

p = inputParser;

% ISOSURFACE VALUE: number at which the levelset of the network function
% will be taken
addOptional(p,'IsoVal',0);

% MISOSURFACE VALUE: number at which the levelset of the network function
% will be taken
addOptional(p,'mIsoVal',4);

% EXTENTS: bounding box for the finger (of the form [-x,x;-y,y;-z,z])
addOptional(p,'Extents',[-4,4;-4,4;0,2*pi]);

% DX: Cell Spacing for Isosurface Evaluation
addOptional(p,'dx',0.05);

% ELEMENTSIZE: determines distance between points used to create the mesh
addOptional(p,'ElementSize',4);

% MEMBRANE THICKNESS: thickness of the shell elements (in mm)
addOptional(p,'MembraneThickness',1.6);

% MEMBRANE INTEGRATION POINTS: thickness of the shell elements (in mm)
addOptional(p,'MemIntPts',5);

% DENSITY: Material density, used for comupting mass (in g/mm^3)
addOptional(p,'Density',2.6e-4);

% PRESSURE: magnitude of the pressure being applied to the inside surface
% of the mesh during simulation
addOptional(p,'P',200*0.00689476);

% BF_OFFSET: distance in mm between the force measurement block and the
% underside of the actuator
addOptional(p,'bf_offset',5);

% REMOVE DIRECTORY: bool determining whether or not to remove the directory
% of the simulations after execution
addOptional(p,'RemoveDir',true);

% MIN STEP SIZE: minimum incremental load step. Use to control sim time.
addOptional(p,'MinStepSize',5e-4);

% CORES: number of cores to run the job on
addOptional(p,'Cores',12)

% REMESH: should we remesh?
addOptional(p,'Remesh',true);

% CLOSEHOLES: should we remesh?
addOptional(p,'CloseHoles',true);

% MAXFACES: toss the mesh if there are more faces than this
addOptional(p,'MaxFaces',1e5);

% KEEP: never throw the mesh away for any reason
addOptional(p,'Keep',false);

% RMISLANDS: should we check for islands?
addOptional(p,'rmIslands',true)

% CAPFACES: should we cap the scalar field faces with large number?
addOptional(p,'CapFaces',0)

parse(p,varargin{:});
Params = p.Results;

end