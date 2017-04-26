

classdef medimage <handle

    %% Medical image class 
    %
    % A class to load medical images and store related data.
    
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

    properties
        
        Label ='';  %% Image label
        header = [];  %% Header data
        vox2mm = [];  %% Vox to mm affine trans. matrix
        vox2tal = [];  %% Vox to talairach or mni transfm. matrix.
        path='';  %% Path to file
        file='';  %% File name
        subjectID='';  %% subject/patient ID
        userdata = [];  %% Other data from user
        notes = '';
        std_orientation = false;
        vox2orig; %% Any transformations on the volume after loading from file.
    end
    
    properties (Dependent=true)
     Data
     load2mem %% Load image data to memory if true, otherwise read from disk.
     vox2std
    end
    
    properties (Hidden=true)
        
        imgDat = [];
        inmem = true;
    end
    
    methods
        
        function me = medimage(varargin)
            
            global MEDIM_LASTPATH
            
            if ~isempty( MEDIM_LASTPATH)
                me.path = MEDIM_LASTPATH;
            end
            
            if nargin == 1 && isempty(varargin{1}) %%% return a volumeview without loading data if the only input is an empty matrix.
                me = me([]);
                return 
            end
            
            if nargin > 0 && isa(varargin{1},'medimage')
                pr = setdiff(properties(varargin{1}),{'vox2std'});
                for k = 1:length(pr)
                    me.(pr{k}) = varargin{1}.(pr{k});
                end
                return
            elseif nargin > 0 && exist(varargin{1},'file')
                [pth,fn,ext] = fileparts(varargin{1});
                me.file = [fn,ext];
                me.path = pth;
            elseif nargin>0
                me.Label = varargin{1};
            end
            
            if nargin > 1
                me.file = varargin{2};
            end
            if nargin > 2
                me.path = varargin{3};
            end
            if nargin <2 && ~exist(fullfile(me.path,me.file),'file')
                me.asgnData;
            else
                me.asgnData(fullfile(me.path,me.file));
            end
             
            me.vox2orig = transforms('trmat',eye(4));
           
            
            MEDIM_LASTPATH = me.path;
            
        end
        
        
        %%%%%%%%
        function asgnData(me,a)
            
            if nargin > 1 && (isnumeric(a) || islogical(a))
                me.imgDat = a;
            else
                if nargin < 2  || isempty(a) || exist(a,'dir')
                     if nargin >= 2
                            me.path = a;
                     end
                     [fn,pth] = uigetfile({'*.nii*','*.hdr';'nifti','img'}',sprintf('Find a %s image for %s',me.Label, me.subjectID),me.path);

                else
                    [pth,fn,ext]= fileparts(a);
                    fn = [fn,ext];
                    if nargin > 2 
                        me.path = pth;
                    end
                end                
                me.path = pth;
                me.file = fn;
                me.loadfile;
            end
        end
            
        function set.Data(me,a)
            asgnData(me,a);
        end
        
        %%%%%%%%%%% 
        function x = get.Data(me)
            
            if me.inmem
                x=me.imgDat;
            else
                [x,me.header] = loadfile(me);
            end
        end
        %%%%%%%
        function set.load2mem(me,a)
            if ~islogical(a) && ~ismember(a,[0 1])
                error('Input must be true (1) or false (0)')
            end
            if isempty(me.imgDat)
                me.inmem = a;
                if a
                    [me.imgDat,me.header] = me.loadfile;
                else
                    me.imgDat=[];
                end
            end
        end
        
        %%%%%%%
        function a=get.load2mem(me)
            a= me.inmem;
        end
        %%%%%%%%%%% load image
        function [img,hdr] = loadfile(me)
              [p,f,ext ]= fileparts(me.file);
              fn = fullfile(me.path,me.file);
               switch lower(ext)
                    case {'.hdr'}
                        [img,hdr] = loadimg(fn);

                    case {'.nii','.gz'}
                          [img,hdr] = readnifti(fn);    
%                    case {'.gz'}
%                        gunzip(fn)
%                        me.file =fullfile(p,f);
%                        me.loadfile();
%                        return
                    otherwise
                    error('Unrecognized file extension %s',ext)
               end
               nd = ndims(img);
               perm = 1:nd;
               e = eye(nd+1); e = e([perm,nd+1],:);
              img = permute(img,perm);
               me.header = hdr;
               me.vox2mm = e(1:3,1:3)*hdr.vox2unit(perm(1:3),:)*e(1:4,1:4)'^-1;
               if me.load2mem
                   me.imgDat = img;
               end
               if length(me.vox2mm)>=3
                   me.std_orientation = all(diag(me.vox2mm(1:3,1:3))==[-1 1 1 ]');
               end           
        end   
        %%%%%
        function tr = get.vox2std(me)
            %Return the standard orientaiton matrix
            stdorient = diag([-1 1 1]);
            curr_orient = sign(me.vox2mm(1:3,1:3)).*(abs(me.vox2mm(1:3,1:3))==repmat(max(abs(me.vox2mm(1:3,1:3))),3,1));
            chorient = curr_orient'*stdorient;
            sz = size(me.Data);
            zer = (sz(1:3)-1)'.*any(chorient<0,2);
%             perm = double(abs(chorient)*(1:3)');
%             if any(zer~=0)
%                 Q = me.Data;
%                 for i = find(zer~=0)
%                    Q=flipdim(Q,i);
%                 end
%                 me.imgDat = Q;
%             end
%             me.imgDat = permute(me.imgDat,perm);
            trmat = chorient;
            trmat(:,4) = zer;
            trmat(4,4) = 1; 
            tr = transforms('trmat',trmat','label','vol2std');
        end
        %%%%%
        function reorient2std(me)
            
            % Reorient the volume into standard orientation
            trmat = me.vox2std.trmat;
            sz = size(me.Data);
            zer = (sz-1).*any(trmat(1:3,1:3)<0);
            %             perm = double(abs(chorient)*(1:3)');
            perm = double(abs(trmat(1:3,1:3))*(1:3)');
            if any(zer~=0)
                Q = me.Data;
                for i = find(zer~=0)'
                   Q=flipdim(Q,i);
                end
                me.imgDat = Q;
            end
            me.imgDat = permute(me.imgDat,perm);
            T = me.vox2mm';
            T(4,4) = 1;
            me.vox2mm = (trmat*T(:,1:3))';
            me.vox2orig = me.vox2orig*transforms('trmat',trmat^-1,'label','std2orig');
            
            
        end
        %%%%%
        function save(me,uigui)
            if nargin <2
                uigui = true;
            end
            if uigui  
                [fn,pth]= uiputfile({'*.nii;*.nii.gz'},'Save file as',fullfile(me.path,me.file));
            else
               fn = me.file;
               pth = me.path;
            end
            
            T = me.vox2mm';
            h = me.header;
            h.srow_x = T(:,1);
            h.srow_y = T(:,2);
            h.srow_z = T(:,3);
            writenifti(h,me.Data,fullfile(pth,fn));
            
        end
            
    end
end
