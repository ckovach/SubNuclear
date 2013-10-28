% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

classdef meshes < plotobj
    properties
            trirep = [];
            file;
            path;
    end
    
     methods
        function imo = meshes(varargin)
            for i = 1:2:length(varargin)
                imo.(varargin{i})=varargin{i+1};
            end
        end
    end
end
