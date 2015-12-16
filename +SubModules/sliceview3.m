
classdef sliceview3 < SubModules.module

%
% The sliceview3 module generates 3 orthogonal planar views of arbitrary orientation. 
%
% To use sliceview3, select 3 or more points or 1 or more meshes containing
% at least 3 vertices in the respective menus, then double click on the 
% 'sliceview' module. A new figure window with three axes will appear.
%
% The orientation of the first view is determined from the first two principle components 
% of the group of selected points, the second by the first and third, and
% the third, by the second and third.
%
%
%
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/SubNuclear/+SubModules/sliceview.m $
% $Revision: 396 $
% $Date: 2013-10-28 10:26:37 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2013

    methods

        function me = sliceview3(vv,X,addsis,varargin)        
            me.initialize(vv,varargin{:});
            vv = me.parent;
            cro = vv.current_object;
            fld = class(cro);
            indx = ismember([vv.(fld)],cro);

            if nargin < 2 || isempty(X)
                
                switch class(cro)
                    case {'points','lines'}


                        X = cat(1,vv.points(indx).coord);

                    case 'meshes'
                        decell = @(x)cat(1,x{:});
                        X = decell( arrayfun(@(x) x.trirep.X, vv.meshes(indx),'uniformoutput',0));
%                     case {'volumes','transforms'}
                        
                    otherwise
                        warning('To use %s first select points or meshes used to define the plane of view.',mfilename)
%                         warning('Showslice does nothing for %s',class(cro))
                        return          

                end
            end
%                 
            if nargin < 3 
                addsis = true;
            end

            if size(X,1) < 3
                error('Requires at least 3 points');
            end

            pt = mean(X);
            vv.current_point = pt;
            [u,~]= svd(cov(X));
%             nv = u(:,3)';
            u(4,4) = 1;

             axsz = [.35 .4];
             axpos = [.05 .55 axsz
                     .55 .55 axsz
                     .3 .05 axsz];
           nfig = figure;
            vv.addfig(nfig);
            cols = 'rby';
            for i = 1:3
                ax(i) = vv.addaxis(axes('position',axpos(i,:),'units','normalized','parent',nfig),false); %#ok<*AGROW>
                ax(i).Transform = transforms('trmat', u(:,[mod((1:3)-i,3)+1, 4]));
                ax(i).setvolmat(ax(i).parent.current_point,false);
                ax(i).crosscol = cols(i);

              
            end
            if addsis
                for k = 1:length(vv.sisters)
                    sis = SubModules.sliceview3(vv.sisters(k),X,false);
                end
                for i = 1:3
                   sis.data.ax(i).sisters = ax(i);
                    ax(i).sisters = sis.data.ax(i);
                   
                end
            end
                vv.resetaxis(ax(:)); 
            
            

           
            me.data.ax = ax;
%             vv.plotupdate();
        end
        
        function update(me)
           % Updating already handled 
        end
            
        
    end
end