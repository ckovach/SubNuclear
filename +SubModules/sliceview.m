
classdef sliceview < SubModules.module

%
% The slice module generates a planar view of arbitrary orientation. 
%
% To use sliceview, select 3 or more points or 1 or more meshes containing
% at least 3 vertices in the respective menus, then double click on the 
% 'sliceview' module. A new figure window and axis will appear.
%
% The orientation of view is determined from the first two principle components 
% of the group of selected points. 
%
% An extra line will appear in the crosshars for all other open axes,
% indicating the lines of intersection for the plane.
%
% You can create a view of arbitrary orientation by selecting 3 points in 
% the desired plane. 
%
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2013

    methods

        function me = sliceview(vv,X,varargin)        
            me.initialize(vv,varargin{:});
            vv = me.parent;
            cro = vv.current_object;
            fld = class(cro);
            indx = ismember([vv.(fld)],cro);

            if nargin < 2 || isempty(X)
                
                switch class(cro)
                    case {'points','lines'}

    
                   X = vv.transforms(1).tr(cat(1,vv.points(indx).coord));


                    case 'meshes'
                        decell = @(x)cat(1,x{:});
                        X = vv.transforms(1).tr(decell( arrayfun(@(x) x.trirep.X, vv.meshes(indx),'uniformoutput',0)));

%                     case {'volumes','transforms'}
                        
                    otherwise
                        warning('To use %s first select points or meshes used to define the plane of view.',mfilename)
%                         warning('Showslice does nothing for %s',class(cro))
                        return          

                end
            end
%                 
%             if nargin < 3 
%                 addsis = true;
%             end

            if size(X,1) < 3
                error('Requires at least 3 points');
            end

            pt = mean(X);
%             vv.current_point = pt;
            vv.current_point = vv.transforms(1).itr(pt);
            [u,~]= svd(cov(X));
            nv = u(:,3)';
            u(4,4) = 1;
%             vc = size(vv.current.volumes.image.Data)/2;
% %              p0 = ((pt-vc)*nv')*nv+vc;
%              sz = size(vv.Vol);
%              %Projection matrix into the plane of view
%             prjmat = [eye(3,4);-pt,1]*(eye(4) - [nv,0]'*[nv,0])*[eye(3,4);+pt,1];
%             
%             %Find the range of values to plot by projecting the corners
%             corn = (crossn([2 2 2])-1)*diag(sz);
%             prjcrn = [corn,ones(size(corn,1),1)]*prjmat;
%             prlim = minmax(prjcrn*[eye(3);-pt]*u(:,1:2));
%             dprlim = diff(prlim);
%             
%             ax2vol = [ [eye(2);prlim(1,:)]*u(:,1:2)',[0 0 1]']*[eye(3);pt];
%             ax2vol(3,4) = 1;
%             vol2ax = ax2vol^-1;
%             axlim = floor(dprlim);
% %             [y,x] = meshgrid(0:1:axlim(2),0:1:dprlim(1));
% %             xy = [x(:),y(:),ones(size(x(:)))];
% %             xytr = [x(:),y(:),ones(size(x(:)))]*ax2vol;
%             
%             
% %             [sl,axv,trmat] = vv.makeslice(nv,p0);
            nfig = figure;
            vv.addfig(nfig);
            axh = axes;
            ax=vv.addaxis(axh); 
            for k = 1:length(ax)
                ax(k).Transform = vv.transforms(1)*transforms('trmat', u);
%                 ax(k).normvec=nv;
%                 ax(k).setvolmat;
                ax(k).setvolmat(vv.transforms(1).tr(ax(k).parent.current_point),false);

            end        
%             szsl = size(sl);
%             axlim(end+1) = 0;
%             ax.axdim = axlim([3 2 3 1]);
            
            % axi = find(ax==vv.axes);

            
            vv.resetaxis(ax); 
%             if addsis
%                 for i = 1:length(vv.sisters)
%                     SubModules.sliceview(vv.sisters(i),X,false);
%                 end
%             end
            me.data.ax = ax;
%             vv.plotupdate();
        end
        
        function update(me)
           % Updating already handled 
        end
            
        
    end
end