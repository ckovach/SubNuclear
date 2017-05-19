
function vtks = readvtk(fns)

%
% Reads vtk files
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2013

persistent pth0
if nargin < 1 || isempty(fns)
    if isnumeric(pth0)
        pth0 ='';
    end
    
   [fname,pth]=uigetfile({'*.vtk;*.bmv;','VTK file';'*.*;','All files'},'Load VTK file',pth0,'multiselect','on');
   if pth0==0
       return
   end
   fns = fullfile(pth,fname);
   pth0=pth;
end

if ~iscell(fns)
    fns = {fns};
end

for kk = 1:length(fns)
    fn = fns{kk};
    fid = fopen(fn,'r');
    vtk = struct('header',{{}});
    for i = 1:4,
        vtk.header{i,1} = fgetl(fid);
    end

    isbinary = isequal(lower(deblank(vtk.header{3})),'binary');


    currfield = 'xxx';
    doread = true;
    while doread
        doread = ~feof(fid);

        l = fgetl(fid);

        if  isequal(l,-1)
            break
        end

        if ~isbinary

            scl = sscanf(l,'%f');

            if isempty(scl)
                scl = sscanf(l,'%s',1);
                switch scl
                    case 'POINTS'
                        currfield = 'vert';
                        vtk.npoints = sscanf(l,'%*s %i %*s');
                        vtk.fmt  = char(sscanf(l,'%*s %*i %s')');
                        vtk.(currfield)(vtk.npoints,1:3) = 0;
                        fun=@(x)x;
                        k = 1;
                    case 'POLYGONS'
                        vtk.n = sscanf(l,'%*s %i %i')';
                        currfield = 'tri';
                        vtk.(currfield)(vtk.n(1)-1,1:3) = 0;
                        fun = @(x)x(2:end)+1;
                        k = 1;
                    otherwise
                        currfield = scl;
                end
    %                 vtk.(currfield) = [];
            else
                vtk.(currfield)(k,:) =fun(scl);
                k = k+1;
            end
        else
                scl = sscanf(l,'%s',1);

                switch scl
                    case 'POINTS'
                        currfield = 'vert';
                        vtk.npoints = sscanf(l,'%*s %i %*s');
                        fmt  = char(sscanf(l,'%*s %*i %s')');
                        fun = @(x)reshape(x,3,length(x)/3)';
                        readn = vtk.npoints*3;
                    case 'POLYGONS'
                        n = sscanf(l,'%*s %i %i')';
                        vtk.ntess = n(1);
                        currfield = 'tri';
                        fmt = 'int32';
                        readn = n(2);
                        rmcol = @(x)x(:,2:end);
                        fun = @(x) rmcol(reshape(x,4,length(x)/4)')+1;
                    case 'FIELD'
    %                     currfield = 'field';
                        fldn = sscanf(l,'%*s %*s %i')';
                        fldi = 1;
                        currfield = char(sscanf(l,'%*s %s %*i')');
                        readn = 0;
                        fmt = 'int';
                        fun=@(x)[];
                    otherwise
                        if fldi <=fldn %#ok<BDSCI>
                            n2 = sscanf(l,'%*s  %i %i %*s')';
                            fnm = sscanf(l,'%s',1); 
                            fmt = char(sscanf(l,'%*s  %*i %*i %s')');
                            readn = prod(n2);
                            fun = @(x) cat(2,vtk.(currfield),struct('label',fnm,'data',reshape(x,n2)));
                            fldi = fldi+1;
                        else
                            error('unrecognized data type')
                        end
                end

                 vtk.(currfield) = fun(fread(fid,readn,fmt));
        end
        
    end
    vtks(kk)=vtk;
    fclose(fid);
end

