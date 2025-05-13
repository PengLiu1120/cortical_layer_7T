function write_vtk(filename,header1,header2,header3,vertex,face,data)

    % Writes data into *.vtk file
    %
    % The function needs outputs from function read_vtk(): header1,
    % header2, header3 and vertex, face, and scalar (or other name)

    fid = fopen(filename,'wt');

    % Write first header information
    for i = 1:numel(header1)
        fprintf(fid,'%s\n', header1{i});
    end
    
    % Write vertex data
    for i = 1:numel(vertex(1,:))
        fprintf(fid, '%.5f %.5f %.5f\n',vertex(:,i));
    end

    % Write second header information
    for i = 1:numel(header2)
        fprintf(fid,'%s\n', header2{i});
    end

    % Write face data
    for i = 1:numel(face(1,:))
        fprintf(fid, '%d %d %d %d\n',face(:,i));
    end

    % Write third header information
    for i = 1:numel(header3)
        fprintf(fid,'%s\n', header3{i});
    end

    % Write scalar data (intensities)
    for i = 1:numel(data)
        fprintf(fid, '%d\n',data(i));
    end

    fclose(fid);
    
return
