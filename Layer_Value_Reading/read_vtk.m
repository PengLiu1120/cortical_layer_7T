function [vertex,face,scalar,header1,header2,header3] = read_vtk(filename)

% Reads data from *.vtk file
%
% 'vertex' is array with all x,y,z coordinates of all vertices  
% 'face' is an array specifying the connectivity of the mesh
% 'scalars' is an vector which specifies the corresponding intensities

fid = fopen(filename,'r');
if( fid == -1 )
    error('Can''t open the file.');
    return;
end

%%% read first header (points) %%%
i = 1;
nvert = 0;
while i > 0
    line = fgetl(fid); 
    header1{i} = line;
    i = i + 1;
    nvert = sscanf(line,'%*s %d %*s', 1);
    if nvert ~= 0
       i = 0;
    end
end

%%% read vertices %%%
[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;

%%% read second header (polygons) %%%
i = 1;
npols = 0;
fgetl(fid);
while i > 0
    line = fgetl(fid);
    header2{i} = line;
    i = i + 1;
    npols = sscanf(line,'%*s %d %*d',1);
    if npols ~=0 
        i = 0;
    end
end

%%% read polygons %%%
[B,cnt] = fscanf(fid,'%d %d %d %d', 4*npols);
if cnt~=4*npols
%     warning('Problem in reading polygons.');
end
B = reshape(B, 4, cnt/4);
face = B;

%%% read third header (scalars) %%%
i = 1;
fgetl(fid);
while i > 0
    line = fgetl(fid);
    header3{i} = line;
    i = i + 1;
    nscal = sscanf(line,'%*s %s',1);
    if strcmp(nscal,'default') == true 
        i = 0;
    end
end

%%% read scalars %%%
[C,cnt] = fscanf(fid,'%f', nvert);
if cnt~=nvert
    warning('Problem in reading scalars.');
end
C = reshape(C, 1, cnt);
scalar = C;

fclose(fid);

return

