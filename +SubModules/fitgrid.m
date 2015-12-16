
classdef fitgrid < SubModules.module

%
% Fitgrid interpolates a grid of specified size using the selected points 
% as control points. 
%
% It accomplishes this using linear affine or non-linear thin-plate spline warping.
%
% The default method is linear if 4 or fewer points are selected (a minum of 3
% is reuired) and non-linear thin-plate spline warping otherwise.
%
% see also FITGRID


% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2013
    properties
        usedefault = true; %Assume defailt grid
        grid_size = [8 12]; %Default grid size
        default_order = [3 3; 1 3; 3 2; 1 2];  %Default point ordering
        transf; 
%         itransf
        points; %Array of point objects
        gridpoints;
        gridscale = [ 5 5]; % Grid scale in mm
    end
    properties (SetAccess = private, Hidden = true)
        fig;
        szh;
        tabh;  
        buth;
        warph; 
        lblh;
        mkaxh;
        rigidspaceh;
    end
    properties (Dependent = true)
%         gridpoints % Points on the grid template (space from which will be warped)
%         size    
        ptlabel
    end
    methods


        function me = fitgrid(varargin)   
            %Fit grid constructor
            me.initialize(varargin{:});
            vv = me.parent;
            pts = vv.current.points;
           
            me.points=pts;

%             lbls = {pts.label};

            
            me.fig = figure('Name','Grid interpolation','numbertitle','off');
            set(me.fig,'units','characters','position',[60   30   50   30])
%             set(me.fig,'DeleteFcn',@(a,b)delete(me))
           
        %                 cpt(1) = uicontrol('style','list','max',2,'units','normalized','position',[.1, .1, .1,.8],'string',lbls);
             sz = uitable('columneditable',true,'units','normalized','rowname','grid size','columnname',{'I','J'});
               pos = [.1, .75 0 0] + get(sz,'extent')*[1 0 0 0; 0 0 0 -1; 0 0 1 0; 0 0 0 1]';
            set(sz,'position',pos,'celleditcallback',@(src,evnt)me.szcallback(src,evnt),'data',me.grid_size)
             me.szh = sz;

              tab = uitable('columneditable',[false,true(1,2)],'units','normalized','data',num2cell(zeros(length(pts),3)),'rowname',{});
              set(tab,'cellselectioncallback',@(src,evnt)me.tabselection(src,evnt))
              pos = [.1, .6 0 0] + get(tab,'extent')*[1 0 0 0; 0 0 0 -1; 0 0 1 0; 0 0 0 1]';
              pos(2) = max(pos(2),.1);
              pos(4) = min(pos(4),.5);
             set(tab,'position',pos,'celleditcallback',@(src,evnt)me.tabcallback(src,evnt),'columnname',{'Point','Grid I','Grid J'})             
             me.tabh = tab;
              
             me.buth = uicontrol('style','pushbutton','units','normalized','position',[.5, .85, .2,.1],'string','OK','callback',@(src,evnt)okcallback(me,src,evnt));
             me.warph(1) = uicontrol('style','radio','units','normalized','position',[.3, .875, .2,.06],'string','Linear','callback',@(src,evnt)wrpcallback(me,src,evnt));
             me.warph(2) = uicontrol('style','radio','units','normalized','position',[.3, .825, .2,.06],'string','TPS','callback',@(src,evnt)wrpcallback(me,src,evnt));
             me.warph(3) = uicontrol('style','radio','units','normalized','position',[.3, .775, .2,.06],'string','Rigid','callback',@(src,evnt)wrpcallback(me,src,evnt));
             me.rigidspaceh = uicontrol('style','edit','units','normalized','position',[.5, .775, .125,.05],'string','[ 5, 5 ]','background',[1 1 1],'callback',@(src,evnt)wrpcallback(me,src,evnt),'enable','off');
             uicontrol('style','text','units','normalized','position',[.1, .9, .2,.05],'string','Label','callback',@(src,evnt)wrpcallback(me,src,evnt));
              me.lblh =uicontrol('style','edit','units','normalized','background',[1 1 1],'position',[.1, .875, .2,.05],'string','Grid %i');
            me.mkaxh  =uicontrol('style','checkbox','units','normalized','position',[.7, .8, .25,.1],'string','Make fig.','value',1);
            
           
             if length(pts)<=4
                set(me.warph(1),'value',1);
                me.wrpcallback(me.warph(1)) 
             else
                set(me.warph(2),'value',1);
                me.wrpcallback(me.warph(2)) 
             end            
             try %%% If any labels of the form x,y are detected, then assume these are grid coordinates.
                ptlbls = {pts.label};
                re = regexp(ptlbls,'(\d+),(\d+)','tokens','once');
                if ~isempty([re{:}]) 
                    re(cellfun(@(x)isempty([x{:}]),re)) = {{'0','0'}};
                    re = cat(1,re{:});
                    grpt = cellfun(@str2double, re);
                    me.gridpoints = grpt;
                    me.usedefault = false;
                    set(me.tabh,'data',cat(2,{me.points.label}',num2cell(grpt)));
      
                end
             catch cerr
                 warning('Caught error ''%s''',cerr.message)
                 me.usedefault = true;
             end
           me.szcallback();

        end
        function tabselection(me,~,evnt)
           if evnt.Indices(2) == 1
              me.parent.current_point = me.points(evnt.Indices(1)).coord;  
           end
        end
%         function a = get.size(me)
%            a =get(me.szh,'data'); 
%         end
%         function a = get.gridpoints(me)
%            X =get(me.tabh,'data'); 
%            a = cell2mat(X(:,2:end));
%         end
         function a = get.ptlabel(me)
           a =get(me.lblh,'string'); 
        end
        function update(me) %#ok<MANU>
            % Nothing to update
        end
        %%%
        function szcallback(me,~,~)
            X = get(me.szh,'data');
            me.grid_size = X;
            if me.usedefault
                X(3) = 1;
                X(4) = 0;
                txmat = me.default_order;
                txmat = txmat(1:min(length(me.points),4),:);
                txmat(min(length(me.points),4)+1:length(me.points),:) = 4;
         
                 set(me.tabh,'data',cat(2,{me.points.label}',num2cell(X(txmat))));
                 me.gridpoints = X(txmat);
            end
        end
        
        function tabcallback(me,~,~)
           me.usedefault = false; 
             X =get(me.tabh,'data'); 
            me.gridpoints = cell2mat(X(:,2:end));
        end
        function wrpcallback(me,src,~) 
%             val = get(me.warph(me.warph==src),'value');
%             set(me.warph(me.warph~=src),'value',1-val);
            set(me.warph(me.warph==src),'value',1);
            set(me.warph(me.warph~=src),'value',0);
            if src == me.warph(3)
               set(me.rigidspaceh,'enable','on') 
            else
                   set(me.rigidspaceh,'enable','off') 

            end
        end
        
        function okcallback(me,~,~) 
                       
           pts = me.points;
           Y = cat(1,pts.coord);
           vox2mm = me.parent.transforms(1);
           Y = vox2mm.tr(Y);
%            X = get(me.tabh,'data');
           warptype = find(cell2mat(get(me.warph,'value')));
             X = me.gridpoints;
           if warptype==3
               gridsc = str2num(get(me.rigidspaceh,'string')); %#ok<ST2NM>
               X = X*diag(gridsc);
           else
               gridsc = me.gridscale;
           end
             [pts,me.transf] = fitgrid(Y,me.grid_size,X,warptype,me.ptlabel,gridsc);     
             pts = vox2mm.itr(cat(1,pts.coord));
           me.parent.addpoint(me.ptlabel,pts);
           
           if get(me.mkaxh,'value')
                
            
                grp = me.gridpoints;
                grsz = me.grid_size ;
                pldim = 5*kron(max(me.grid_size)*[1 1],[-1 2]);
                [x,y] = meshgrid(pldim(3):pldim(4),pldim(1):pldim(2));
             
                [xgr,ygr] = meshgrid((1:grsz(2))*5,(1:grsz(1))*5);

              
                xyi = vox2mm.itr(me.transf.tr([y(:),x(:)]));

                q =me.parent.current.volumes.interpolant(xyi(:,1),xyi(:,2),xyi(:,3));

                % xyl = sub2ind(size(vv.Vol),xyi(:,1),xyi(:,2),xyi(:,3));


                q = reshape(q,size(x));
                %%%
                figure, imagesc(x([1 end]),y([1 end]),q);
                caxis(me.parent.intensity_range)
                colormap gray
                 colormap gray, axis xy, 
                hold on, plot(xgr(:),ygr(:),'b.')
                hold on, plot(grp(:,2)*5,grp(:,1)*5,'r+')
                axis image
%                axh = axes;
%                ax = me.parent.addaxis(axh); 
%                ax.Transform = itrans;
%                
%                ax.axdim = [[-.25 1.25].* me.grid_size(2), [-.25 1.25].*me.grid_size(1)]*diag(gridsc([2 2 1 1]));
%                ax.setvolmat;
%                ax.parent.plotupdate(ax);
           end
           if ~strcmp(get(me.fig,'beingdeleted'),'on')
               delete(me.fig);
           end
           
        end
    end
end