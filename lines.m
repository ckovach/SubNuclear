% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

classdef lines < points
     methods
        function imo = points(varargin)
            for i = 1:2:length(varargin)
                imo.(varargin{i})=varargin{i+1};
            end
        end
    end
end
