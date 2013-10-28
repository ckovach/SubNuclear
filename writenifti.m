
function writenifti(hdat,IM,fname)

%function writenifti(hdat,IM,fname)
%Write nifti file, return image and header data. Requires the format specification file,
%nfti_hdr_fmt.txt to be present.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2010

fid = fopen(fname,'w','l');

% hdat.endian = 'little';
 fwrite(fid,hdat.sizeof_hdr,'int32');
% if ~isequal(hdat.sizeof_hdr,348) % hdat.sizeof_hdr should always be 348. If that's not the case then we are either using the wrong endianess or the file format is wrong.
%     fclose(fid);
%     fid = fopen(fname,'r','b');  %Try big endian
%     hdat.sizeof_hdr = fread(fid,1,'int32');
%     if ~isequal(hdat.sizeof_hdr,348)
%         error('file format is not as expected')
%         fclose(fid)
%     end
% 
%     hdat.endian = 'big';
% end

datafmts = {'ubit1','uchar','int16','int32','single','double','double[2]','double[3]','char','uint16','uint32','int32','uint64','float128','double[2]','float128[2]'};


hdat.dim1 = size(IM,1);
hdat.dim2 = size(IM,2);
hdat.dim3 = size(IM,3);
hdat.dim4 = size(IM,4);
hdat.dim5 = size(IM,5);
hdat.dim6 = size(IM,6);
hdat.dim7 = size(IM,7);
    
%%
hdr_fmt = ReadTabDelim('nfti_hdr_fmt.txt','\t');
for i = 3: size(hdr_fmt,1)
    
    k = i-1;
    
    typedat = regexp(hdr_fmt{i,4},'(\w*)\[?(\d*)\]?','tokens');
    typedat = [typedat{:}];
    if isempty(typedat{2}), typedat{2} = '1'; end
    fmt(k).label= deblank(hdr_fmt{i,3});
    fmt(k).type = typedat{1};
    fmt(k).number= str2num(typedat{2}); %#ok<*ST2NM>
    fmt(k).bytes= str2num(hdr_fmt{i,2}); %#ok<*AGROW>
    fmt(k).info= hdr_fmt{i,5};
    
    fn = regexprep(fmt(k).label,'[\[\]]','');
%     hdat.(fn) = fread(fid,fmt(k).number,[fmt(k).type,'=>',fmt(k).type]);
    fwrite(fid,hdat.(fn),fmt(k).type);
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
fwrite(fid,[0 0 0 0],'uchar');
stat = fseek(fid,hdat.vox_offset,-1);
if stat <0
   error(ferror(fid)); 
end
fwrite(fid,IM(:),format);
fclose(fid);

[~,~,ext] = fileparts(fname);


switch lower(ext)
    case {'.gz'}  %%% Zip file if it ends with gz
        gzip(fname)
end




%%


  

