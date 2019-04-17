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
            if nargin == 1 && isempty(varargin{1})
                imo = imo([]);
                return
            end
            if nargin> 0 && isa(varargin{1},'meshes')
                fldn = fieldnames(varargin{1});
                fldn = setdiff(fldn,{'ploth','objectid'});
                for k = 1:length(fldn)
                    imo.(fldn{k}) = varargin{1}.(fldn{k});
                end
            elseif nargin>0
                
                for i = 1:2:length(varargin)
                    imo.(varargin{i})=varargin{i+1};
                end
            end
        end
    end
end
