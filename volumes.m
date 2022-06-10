% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

classdef volumes < plotobj
    properties
           intensity_range = [0 255];
           file;
           path;
           interpolant;
           parent;
    end
    
    properties (Dependent = true)
        image;
    end
    properties (Hidden = true, SetAccess = private)
        medimg = medimage([]);
    end
     methods
        function imo = volumes(varargin)
            for i = 1:2:length(varargin)
                imo.(varargin{i})=varargin{i+1};
            end
        end
        function a = get.file(me)
            a = me.image.file;
        end
        function a = get.path(me)
            a = me.image.path;
        end
        
        function a = get.interpolant(me)
            if isempty(me.interpolant) 
                X = me.image.Data;
                me.interpolant = griddedInterpolant({0:size(X,1)-1,0:size(X,2)-1,0:size(X,3)-1},X,'linear','none');
            end
            a = me.interpolant;
        end
        
         function tr = tr2std(me)
            % Returns transform to standard orientation based on sign of
            % vox2mm 
            tr = me.image.vox2std;
                                    
            
         end
         function reorient2std(me)
             
             me.image.reorient2std;
             if isa(me.parent,'volumeview')
                me.parent.transforms(1).trmat = me.image.vox2mm';
                me.parent.transforms(2) = me.tr2std;
             end
             me.interpolant = [];
         end
         function set.image(me,a)
             me.medimg = medimage(a);
             me.interpolant = [];
         end
         function a = get.image(me)
            a = me.medimg;
         end
     end
end
