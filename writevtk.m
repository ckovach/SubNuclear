
function writevtk(vtk,fn)


% Writes mesh data to an ASCII vtk file. 

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2013

fid = fopen(fn,'w');

vtk.header{3,1} = 'ASCII';
vtk.header{2,1} = 'vtk output';
vtk.header{1,1}='# vtk DataFile Version 3.0';
vtk.header{4,1}='DATASET POLYDATA';
vtk.fmt = 'float';
for i = 1:4,
    fprintf(fid,vtk.header{i,1});
    fprintf(fid,'\n');
end

% isbinary = isequal(lower(deblank(vtk.header{3})),'binary');
vtk.npoints = size(vtk.vert,1);
vtk.n = size(vtk.tri,1)*[1 size(vtk.vert,2)+1];

fprintf(fid,'POINTS %i %s\n',vtk.npoints,vtk.fmt);
fprintf(fid,'%0.9f  %0.9f  %0.9f\n',vtk.vert');

fprintf(fid,'POLYGONS %i %i\n',vtk.n);
ndim = size(vtk.tri,2);
str = sprintf('%i %s\n', ndim,repmat('%i ',1,ndim));

tr = vtk.tri';
tr = tr(:);
% tr(end+1) = 0; %Not sure why this is necessary
fprintf(fid,str,tr(:)-1);

fclose(fid);

    