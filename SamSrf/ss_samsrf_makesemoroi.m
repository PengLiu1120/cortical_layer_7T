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
% 07/10/2021 - Adapted from MakeOccRoi.m (SS)
% 08/02/2022 - Last modified (SS)

Y1 = 5; % anterior-posterior
Y2 = -55; % anterior-posterior
Z1 = 30; % superior-inferior
Z2 = 80; % superior-inferior
XLeft = -2; % right-left
XRight = 2; % right-left
Hemis = {'lh' 'rh'};
% Loop thru hemispheres
for h = 1:size(Hemis,2)
    CurrHemi = Hemis{1,h};
    Vs = fs_read_surf([MeshFolder filesep Hemis{h} '.inflated']);
    Srf = struct;
    Srf.Vertices = Vs;
    Srf.Data = zeros(1,size(Srf.Vertices,1));
    Srf.Hemisphere = CurrHemi;
    if strcmp(CurrHemi, 'lh')
        vx = find(Vs(:,2) <= Y1 & Vs(:,2) >= Y2 & Vs(:,3) >=Z1 & Vs(:,3) <=Z2 & Vs(:,1) <= XLeft);
    else
        vx = find(Vs(:,2) <= Y1 & Vs(:,2) >= Y2 & Vs(:,3) >=Z1 & Vs(:,3) <=Z2 & Vs(:,1) >= XRight);
    end
    samsrf_srf2label(Srf, [CurrHemi '_sensorimotor'], 1, vx);
end

