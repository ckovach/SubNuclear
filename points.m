% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

classdef points < plotobj
    properties
            coord = [];
            file;
            path;
            units = 'vox';
    end
     methods
        function imo = points(varargin)
            
            if nargin == 1 && isempty(varargin{1})
                imo = imo([]);
                return
            end
            if nargin > 0 && isnumeric(varargin{1})            
                X = varargin{1};
                varargin(1) = [];
                if size(X,1) == 1
                    imo.coord = X;
                else
                    for i = 1:size(X,1)
                        imo(i) = points(X(i,:),varargin{:});
                    end
                end
             end
             if nargin > 0 && ~isempty(varargin)>0 && isa(varargin{1},'points')
                fln = fieldnames(varargin{1});
                
                for i = 1:length(fln)
                    for ii = 1:length(varargin{1})                     
                        imo(ii).(fln{i}) = varargin{1}(ii).(fln{i}); 
                    end
                end
             else
                for i = 1:2:length(varargin)
                   if isfield(imo,varargin{i})
                        imo.(varargin{i})=varargin{i+1};
                   else
                       warning('%s is not a valid field name for a %s object.',varargin{i},class(imo))
                   end
                end
            end
        end
    end
end
