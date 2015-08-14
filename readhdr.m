
function [hdat,fmt] = readhdr(fname)


fid = fopen(fname,'r');

fdat = dir(fname);

datafmts = {'ubit1','uchar','int16','int32','single','double','double'};
orientations = {'R-L','P-A','I-S'
            'R-L','I-S','P-A'
            'P-A','I-S','R-L'
            'R-L','A-P','I-S'
            'R-L','S-I','P-A'
            'P-A','S-I','R-L'};
    
%%
hdr_fmt = ReadTabDelim('hdr_fmt.txt');
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


for i = 2: size(hdr_fmt,1)
    
    k = i-1;
    
    typedat = regexp(hdr_fmt{i,4},'(\w*)\[?(\d*)\]?','tokens');
    typedat = [typedat{:}];
    if isempty(typedat{2}), typedat{2} = '1'; end
    fmt(k).label= deblank(hdr_fmt{i,3});
    fmt(k).type = typedat{1};
    fmt(k).number= str2num(typedat{2});
    fmt(k).bytes= str2num(hdr_fmt{i,2});
    fmt(k).info= hdr_fmt{i,5};
    
    fn = regexprep(fmt(k).label,'[\[\]]','');
    hdat.(fn) = fread(fid,fmt(k).number,[fmt(k).type,'=>',fmt(k).type]);
    if hdat.sizeof_hdr ~= fdat.bytes %% Assume this is a problem with endianness
        fclose(fid);
        fid = fopen(fname,'r','b');
        hdat.(fn) = fread(fid,fmt(k).number,[fmt(k).type,'=>',fmt(k).type]);
        if  hdat.sizeof_hdr ~= fdat.bytes
            fclose(fid);
            fid = fopen(fname,'r','l');
            hdat.(fn) = fread(fid,fmt(k).number,[fmt(k).type,'=>',fmt(k).type]);
            if  hdat.sizeof_hdr ~= fdat.bytes %% Go back to native encoding if it's still not working
                fclose(fid);
                fid = fopen(fname,'r','n');
                hdat.(fn) = fread(fid,fmt(k).number,[fmt(k).type,'=>',fmt(k).type]);
            end
        end
    end
end
hdat.dim_orientation = orientations(hdat.orient+1,:);
hdat.vox2unit = cat(2,diag([hdat.pixdim1,hdat.pixdim2,hdat.pixdim3]),zeros(3,1));

%%
fclose(fid)

  
%%

