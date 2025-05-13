function ss_samsrf_makesemoroi(MeshFolder)
%
% MakeSensorimotorRoi(MeshFolder)
%
% Creates a sensorimotor ROI for pRF analysis for each hemisphere.
% The input argument is the MeshFolder (e.g. '..\surf'). 
%
% Note that the vertex coordinates saved in the label file are 
%  from the spherical mesh, -not- the white matter mesh as normally.
%
% 05/10/2021 - Adapted from MakeOccRoi.m (SS)
% 24 Jan 22 / Last modified (PL)

% Y1 = +35; % anterior-posterior
Y1 = +15; % anterior-posterior
Y2 = -55; % anterior-posterior
Z = 10; % superior-inferior

Hemis = {'lh' 'rh'};
% Loop thru hemispheres
for h = 1:2
    Vs = fs_read_surf([MeshFolder filesep Hemis{h} '.inflated']);
    Srf = struct;
    Srf.Vertices = Vs;
    Srf.Data = zeros(1,size(Srf.Vertices,1));
    vx = find(Vs(:,2) <= Y1 & Vs(:,2) >= Y2 & Vs(:,3) >=Z);
    samsrf_srf2label(Srf, [Hemis{h} '_sensorimotor'], 1, vx);
end

