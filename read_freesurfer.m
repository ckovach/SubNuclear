function out = read_freesurfer(filename)

% Read freesurfer binary file

%C Kovach 2013


fid = fopen(filename);
if fid < 0
    error(sprintf('Error opening file'))
end

TRIANGLE_FILE_MAGIC_NUMBER  =  16777214;
QUAD_FILE_MAGIC_NUMBER      =  16777215;
QUAD_FILE_MAGIC_NUMBER2      =  16777213;


%%% magic number identifying format
m = fread(fid,1,'ubit24=>double',0,'b');
%%% NOT YET IMPLEMENTED
switch m

    case TRIANGLE_FILE_MAGIC_NUMBER  %% Binary quadrangle
        out.type = 'triangle';
        out.creator = fgetl(fid);
        fread(fid,1,'uchar'); %extra \n
        out.npoints = fread(fid,1,'int32=>double',0,'b');
        out.n(1) = fread(fid,1,'int32=>double',0,'b');
        vert = nan(3,out.npoints);
        vert(:) = fread(fid,out.npoints*3,'float=>double',0,'b');
        out.vert = vert';
        tri= nan(3,out.n(1));
        tri(:) = fread(fid,out.n(1)*3,'int32=>double',0,'b') + 1;
        out.tri = tri';
        therest = fread(fid,'uchar=>char')';
        voxelsize = regexp(therest,'voxelsize\s*=([^\n]*)','tokens','once');
        voxelsize = str2num(voxelsize{1});
        voldim = regexp(therest,'volume\s*=([^\n]*)','tokens','once');
        voldim = str2num(voldim{1});
        trmatx = regexp(therest,'[xyzc]ras\s*=([^\n]*)','tokens');
        trmatx = cellfun(@str2num,[trmatx{:}],'uniformoutput',false);
        trmatx = cat(1,trmatx{:});
%         trmatx(4,1:3) = voldim/2.*sign(trmatx(4,:));
        out.coord2ras = trmatx;
    case {QUAD_FILE_MAGIC_NUMBER,QUAD_FILE_MAGIC_NUMBER2}  %% Ascii format
        out.type = 'quad';
otherwise
        out.type = '';
end



    
