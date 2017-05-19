
function varargout = readnifti(fname, headeronly)


% [img,hdr] = readnifti(fname)
% Read nifti file, return image and header data. Requires the format specification file,
% nfti_hdr_fmt.txt to be present.
%
% hdr = readnifti(fname,true) 
%
% Returns header data only if the second argument is true.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2010

if nargin < 2
    headeronly = false;
end

[~,~,ext] = fileparts(fname);

deletefile = false; % Delete file after reading. You probably want this to be false unless you work for the CIA.
if strcmp(ext,'.gz')  % If gzipped, check if unzipped file exists and unzip if it doesn't
        tmpfname = [tempname,'.nii.gz'];
        copyfile(fname,tmpfname)
        gunzip(tmpfname);
        deletefile = true; %Unzipped files are deleted to avoid redundant copies of the same data.
        fname = tmpfname(1:end-3);
%     else
%         warning('!!! An unzipped file with the same name already exists in the directory: !!!\n\t%s\n This file will be opened instead. Rename or remove the file if this is not desired.',fullfile(pth,fn))
%     fname = fullfile(pth,fn); 
end

    
fid = fopen(fname,'r','l');

hdat.endian = 'little';
hdat.sizeof_hdr = fread(fid,1,'int32');
if ~isequal(hdat.sizeof_hdr,348) % hdat.sizeof_hdr should always be 348. If that's not the case then we are either using the wrong endianess or the file format is wrong.
    fclose(fid);
    fid = fopen(fname,'r','b');  %Try big endian
    hdat.sizeof_hdr = fread(fid,1,'int32');
    if ~isequal(hdat.sizeof_hdr,348)
        error('file format is not as expected')
        fclose(fid) %#ok<UNRCH>
    end

    hdat.endian = 'big';
end

datafmts = {'ubit1','uchar','int16','int32','single','double','double','double','char','uint16','uint32','int32','uint64','float128','double[2]','float128[2]'};

    
%%
hdr_fmt = ReadTabDelim('nfti_hdr_fmt.txt','\t');
for i = 3: size(hdr_fmt,1)
    
    k = i-1;
    
    typedat = regexp(hdr_fmt{i,4},'(\w*)\[?(\d*)\]?','tokens');
    typedat = [typedat{:}];
    if isempty(typedat{2}), typedat{2} = '1'; end
    fmt(k).label= deblank(hdr_fmt{i,3}); %#ok<*AGROW>
    fmt(k).type = typedat{1};
    fmt(k).number= str2num(typedat{2}); %#ok<*ST2NM>
    fmt(k).bytes= str2num(hdr_fmt{i,2});
    fmt(k).info= hdr_fmt{i,5};
    
    fn = regexprep(fmt(k).label,'[\[\]]','');
    hdat.(fn) = fread(fid,fmt(k).number,[fmt(k).type,'=>',fmt(k).type]);
end
% hdat.dim_orientation = orientations(hdat.orient+1,:);

%% Get transformation into scanner coordinates...

if hdat.qform_code == 0  % Method 1 (ANALYZE compatible)
    T = eye(4);
elseif hdat.sform_code == 0; % Method 2
    offset = [hdat.qoffset_x,hdat.qoffset_y,hdat.qoffset_z]';
    pxscale= [hdat.pixdim1,hdat.pixdim2,hdat.pixdim3]';
    
    %Convert quaternions to rotation matrix....
    Q2R = @(b,c,d)         [ 1-(b^2+c^2+d^2)+b^2-c^2-d^2   2*b*c-2*sqrt(1-(b^2+c^2+d^2))*d       2*b*d+2*sqrt(1-(b^2+c^2+d^2))*c     ;...
                             2*b*c+2*sqrt(1-(b^2+c^2+d^2))*d       1-(b^2+c^2+d^2)+c^2-b^2-d^2   2*c*d-2*sqrt(1-(b^2+c^2+d^2))*b     ;...
                             2*b*d-2*sqrt(1-(b^2+c^2+d^2))*c       2*c*d+2*sqrt(1-(b^2+c^2+d^2))*b       1-(b^2+c^2+d^2)+d^2-c^2-b^2 ];

    R = Q2R(hdat.quatern_b,hdat.quatern_c,hdat.quatern_d)*diag([1 1 hdat.qfac]);
                             
    T = cat(2, R*diag(pxscale) , offset);
    
else
    
    T = cat(1,hdat.srow_x',hdat.srow_y',hdat.srow_z');
    
end

hdat.vox2unit = T;


if hdat.datatype == 0
   warning('Unknown data type')
   lbl = strcat(cellfun(@num2str,num2cell(1:length(datafmts)),'uniformoutput',false),{'. '},datafmts);
   inp = input(['\nAssume which type?',sprintf('\n%s',lbl{:}),'\n']);
   format = datafmts{inp};
else
    
    format = datafmts{log(double(hdat.datatype))./log(2)+1};
end


% if double(hdat.extentder(1))>0
%     %Add code to get extended data here...
% end

if ~headeronly
    fseek(fid,hdat.vox_offset,-1);
    img = fread(fid,format);
    imgdim = double([hdat.dim1, hdat.dim2, hdat.dim3, hdat.dim4, hdat.dim5, hdat.dim6, hdat.dim7] );

    IM = reshape(double(img),imgdim(imgdim>0));
    varargout{1} = IM;
    varargout{2} = hdat;
    varargout{3} = fmt;
    
else
   varargout{1} = hdat;
   varargout{2} = fmt;
end

fclose(fid);

if deletefile
    delete(fname)
    delete(tmpfname)
end

%%


  

