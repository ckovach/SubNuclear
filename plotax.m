
classdef plotax < handle
% Class definition of the plotax object for volume view. It contains
% information about the plane of view for a given axis.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2013
    properties
        label = '';
        notes= '';
        AXV;
        h = -1;
        objectid = -1;
        handles = struct('control',[],'transpose',[],'rotate',[],'crosshairs',[]);
%         NonLinearProjection =false; % True if projection is non-affine
%         normvec;
%         axvec;
        ploth; %  images handles
        rot = 0;
        transpose = 0;
        parent;
        axdim;
%         sisters;
        Transform = transforms('trmat',eye(4));  % Transform
        showcrossHere = true; %Show crosshairs on this axis
        showcrossThere = true; %Show crosshairs on other axes (if showcrossHere is true for them)
        crossprop = '-'; %Line properties of view intersection as plotted in the crosshairs of other axes
        crosscol = 'g';  %Line color of view intersection as plotted in the crosshairs of other axes
        scalebar = struct('length',10,'color',[1 1 1],'position',[.1 .05],'width',.25,'handle',[]);
%         a2vmat;
%         v2amat;
        axtransform = transforms('trmat',eye(4)); % Any additional affine transformation following the transformation into the 2 view plane
        plottol=4; % how many voxels a point can be off the plain and still rendered
    end
     properties (Hidden = true)
         basetr = eye(4);
         normv;
         sisobj;
%          gridpts;
     end
     properties (Dependent = true)
         basetrmat;
         normvec;
         gridpoints;
         trmat;
         sisters;
     end
     
   methods
        
        function me = plotax(varargin)
            
            % me = plotax(parent_volumeview_object, axis_handle)
            %
            %
            
            if nargin > 0 &&  isa(varargin{1},'plotax')
                
                props = setdiff(properties(varargin{1}),'gridpoints');
                for i = 1:length(props)
                   me.(props{i}) = varargin{1}.(props{i}); 
                end
                
            elseif nargin > 0 && ishandle(varargin{1})

                me.h = varargin{1};
            elseif nargin > 0 && isa(varargin{1},'volumeview')
                me.parent = varargin{1};
                if nargin > 1
                    me.init(varargin{2}); 
                else
                    
                    me.init(axes('parent',figure));
                end
%                 me.parent.plotax(end+1) = me;
            end
            
            
        end
        %%%
        
        function set.basetrmat(me,a)
%             
%             me.Transform = transforms('trmat',a);
%             if ~isempty(me.sisters)
%                 me.sisters = me.sisters(isvalid(me.sisters));
%                 for i = 1:length(me.sisters)
%                     me.sisters(i).Transforms= me.Transform;
%                 end
%             end
            me.basetr = a;
            
        end
        %%%
        function cleanup(me)
           me.sisters = unique(me.sisters(isvalid(me.sisters)));
           
        end
        %%%
        function b=get.basetrmat(me)
%             b = me.basetr;
            b = me.basetr;
        end
            %%%
        
%         function set.normvec(me,a)
% %             me.normv = a;
% %             for i = 1:length(me.sisters)
% %                 me.sisters(i).normv= a;
% %             end
%         end
        %%%
        function b=get.normvec(me)
%             b = me.normv;
              b = me.Transform.trmat(1:3,3)';
              
        end 
        %%%
        function ish = isopen(me)
            ish =isvalid(me);
            ish(ish) = ishandle([me(ish).h]);
        end
        %%%
        function c = children(me)
           c = get(me.h,'children'); 
        end
        function out= get(me,varargin)
            out = get(me.h,varargin{:});
        end
        function set(me,varargin)
            set(me.h,varargin{:});
        end
      
        function  setvolmat(ax,pt,updatesis)
            if nargin < 2
               pt = ax.parent.current_point; 
            end
            if nargin < 3
                updatesis = true;
            end
            if any(~isnumeric(pt))  || any(isnan(pt))||any(isinf(pt))
                warning('Non-numeric value detected for current point')
                return
            end
%             T = ax.trmat;
            T = ax.Transform.trmat;
            T(4,1:3) = 0;
            T(4,4) = 1;
            prjmat = [eye(3,4);[-pt,1]]*T; % Projection into plane of view           
            prjmat(1:end-1,end)=0;
            prjmat(end,end)=1;
            
            sz = size(ax.parent.Vol);
            %Find the range of values to plot by projecting the corners
            corn = (crossn([2 2 2])-1)*diag(sz-1);
            corn(:,4) = 1;
%             corn(:,4) = 1;

%             cornax = corn*[eye(3,4);[-pt,1]]*U;
            cornax = corn*prjmat;
%             cornax = ax.vol2ax(corn);
            prlimax = minmax(cornax);
            rg = prlimax(:,3);
            prlimax(:,3) = 0;
            prlimax(:,4) = -1;
           
            vol2ax = prjmat*[eye(3,4);-prlimax(1,:)];
%             vol2ax = prjmat;
            if any(isnan(vol2ax(:)))
                return
            end
%             ax.Transform = transforms('trmat',vol2ax);
% %             ax.normvec = vol2ax(1:3,3);
            ax.Transform.trmat = vol2ax*ax.axtransform.trmat;
            ax2vol = vol2ax^-1;
            
%              ax.basetrmat = vol2ax;
%              ax.v2amat = vol2ax;
             
     
            
            if any(isnan(ax2vol(:)))
                warning('Nan value detected in axis project matrix. Not updating.')
                return
            end
            axlim = floor(diff(prlimax(:,1:2)));
            axlim = min(abs(axlim),ones(size(axlim))*4e3).*sign(axlim);
            axlim(end+1) = 1;
            axlim = ax.axtransform.tr([zeros(1,3);axlim]);
            
            set(ax.handles.control,'Min',rg(1),'Max',rg(2));
            set(ax.handles.control,'Value',0);
            
%             ax.axdim = axlim([3 1 3 2])-1;
            ax.axdim = round(axlim([1 2 3 4])-1);
        
            if get(ax.parent.fixSisterAx,'value') && updatesis
                if ~isempty(ax.sisters)
                    ax.sisters = ax.sisters(isvalid(ax.sisters));
                end
                
                %Make sure sister axes show the same view 
%                 if ~isempty(ax.sisters)
%                     ax.sisters = ax.sisters(isa(ax.sisters,'plotax'));
%                     ax.sisters = ax.sisters(isvalid(ax.sisters));
%                 end
                if ~isempty(ax.sisters)
                    for i = find(isvalid(ax.sisters))

                       ax.sisters(i).Transform = ax.Transform;
                       ax.sisters(i).rot = ax.rot;
                       ax.sisters(i).transpose = ax.transpose;
                       ax.sisters(i).setvolmat(pt,false);
                    end

                end
            end
            
        end
    %%%
     function ptax  = vol2ax(ax,ptvol)
       % Project volume coordinate into axis space
%             vol2ax = ax.trmat;
          ptax = ax.Transform.tr(ptvol);
          ptax(:,4) = 1;
          ptax = ptax*ax.basetrmat(:,1:2);
%           ptax = ptax(:,1:2);
%           T = ax.trmat;
%           ptax = [ptvol,ones(size(ptvol,1),1)]*T(:,1:2);
     end
      %%%
     function ptvol = ax2vol(ax,ptax)
        
         % Project axis coordinate into volume space
         ptax(:,3) = 0;
         ptax(:,4) = 1;
         ptvol = ax.Transform.itr(ptax*ax.basetrmat^-1);
%         T = ax.trmat^-1;
%         ptvol = [ptax,zeros(size(ptax,1),1),ones(size(ptax,1),1)]*T(:,1:3);
        
         
     end
     
     function xytr = get.gridpoints(me)

        if isempty(me.axdim)
            xytr = [];    
        else
            vv = me.parent;
            me.setvolmat(vv.current_point);
            [y,x] = meshgrid(me.axdim(3):me.axdim(4),me.axdim(1):me.axdim(2));
            xytr = me.ax2vol([x(:),y(:)]);
            
        end
         
     end
         
     function trbutton(ax,src)
            
            
            switch src
                case ax.handles.transpose % transpose
                    ax.transpose = get(src,'value'); 
                    vc = ax.vol2ax(size(ax.parent.Vol)/2-.5);
                    ax.basetrmat = ax.basetrmat*[eye(3,4);[-vc,0,1]]*[0  1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]*[eye(3,4);[+vc,0,1]];
                case ax.handles.rotate % Rotate 90 degrees
                    ax.rot=mod(ax.rot+1,4);
                    vc = ax.vol2ax(size(ax.parent.Vol)/2-.5);
                    ax.basetrmat = ax.basetrmat*[eye(3,4);[-vc,0,1]]*[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]*[eye(3,4);[+vc,0,1]];
            end
              ax.setvolmat;
              ax.parent.resetaxis(ax);
           
     end
     function T=get.trmat(ax)
            %Return the rotated transformation matrix             
            T = ax.Transform.trmat*ax.basetrmat;            
     end
     function set.trmat(ax,a)
            ax.Transform = transforms('trmat',a);
%             ax.basetrmat=a;            
     end
     function update(me)
         me.parent.plotupdate(me);
     end
     function delete(me) 
         me.rmaxis;
         if ~isempty(me.parent)
            me.parent.cleanup;
         end
     end
     
     function ax = init(ax,axh)
        
        
        if nargin < 2
            axh = ax.h;
        else
            ax.h = axh;
        end
        vv = ax.parent;
        hold(axh,'on');              

        set(axh,'units','normalized')
        axis(axh,'xy')
        axis(axh,'off')
        axpos = get(axh,'position');
        u =uicontrol( 'style','slider', 'units','normalized', 'position',[axpos(1:2)+[axpos(3)+.005, 0], .03, axpos(4)-.05],'Callback',@ ax.sliderupdate);%,...
%                          'Min',0,'Max',max(size(me.Vol))-1,'SliderStep',[1/max(size(me.Vol)) .1]); 
         utrs = uicontrol( 'style','togglebutton', 'units','normalized', 'position',[axpos(1:2)+[axpos(3)-.0175 axpos(4)-.05], .03, .05],...
                                        'string','T','fontsize',12,'Callback',@(a,b)ax.trbutton(a),'value', false);
         urt = uicontrol( 'style','pushbutton', 'units','normalized', 'position',[axpos(1:2)+[axpos(3)+.0175 axpos(4)-.05], .03, .05],...
                                        'string','R','fontsize',12,'Callback',@(a,b)ax.trbutton(a),'value', false);
       setappdata(utrs,'axis',axh)
       setappdata(urt,'axis',axh)
       
%         me.sliders = [me.sliders,u];
%         ax.normvec=[0 0 1];
%         ax.axvec=[1 0 0; 0 1 0];
        ax.handles.control=u;
        ax.rot=0;
%         ax.handles.trbuttons=[utrs urt];
        ax.handles.transpose=utrs;
        ax.handles.rotate = urt;
%         plh = plot(axh,nan(2,(length(ax.parent.plotax)-1)),nan(2,(length(ax.parent.plotax)-1)),'color','g');
%         ax.handles.crosshairs=plh;
        set(axh,'ButtonDownFcn',@(a,b)ax.axisupdate(a,b),'DeleteFcn',@(a,b)ax.rmaxis(a,b));
        nax=sum( ishandle([vv.plotax.h]) & ax~=vv.plotax);
        ax.handles.crosshairs = plot(ax.h,nan(2,nax),nan(2,nax),'color','g');
        
        for k = find( ishandle([vv.plotax.h]) & ax~=vv.plotax) 
%             xkk = me.plotax(k).crosshairs;
            vv.plotax(k).handles.crosshairs(end+1) = plot(vv.plotax(k).h,nan(2,1),nan(2,1),'color','g');
%             me.plotax(k).crosshairs =xkk;
        end
%         if  ax.start_transposed
%           
%                ax.trbutton(utrs);
%      
%         end 
%         
        set(axh,'DeleteFcn',@(a,b)delete(ax))
        

     end
      
         %%%
     function sliderupdate(me,src,evnt) %#ok<*INUSD>
        
         set(me.handles.control,'Enable','inactive')
         crp = me.parent.current_point;
           u = me.handles.control;
           nv = me.normvec;
%           crp = crp - (nv*crp')*nv + nv*get(u,'Value'); 
          crp = crp  + nv*get(u,'Value'); 
        me.parent.current_point = crp;
%         pause(.1)
         set(me.handles.control,'Enable','on')
       
     end
    %%%%
     function axisupdate(ax,src,evnt)
        axh =ax.h;
        try
            axcp = get(axh,'CurrentPoint');

            ax.parent.current_point=ax.ax2vol(axcp(1,1:2));
        catch %#ok<*CTCH>
           ax.parent.cleanup; 
        end
     end
     
      %%%%%
    function rmaxis(ax,src,evnt)
       if ~ax.isvalid
           return
       end
        flds = ax.parent.object_types;
        for i = 1:length(flds)
            fld =flds{i};
            for k = 1:length(ax.(fld))
                if isa(ax.(fld),'plotobj')
                    ax.(fld)(k).ploth(ax.plotax== ax.h) = [];
                end
            end
        end
        ax.parent.plotax(ax.parent.plotax == ax) = [];
        ax.parent.cleanup;
        
    end
    
    
    
    %%%%%
    
    function plotupdate(me,updatesis)
          
        if nargin < 2
            updatesis = true;
        end
          % Update the plot on this axis
          
          %%% Parent object
          vv = me.parent;
          
          %%% Active volume
          vol =  me.parent.current.volumes(1);
          
          %%% Plotting points
           grxy = me.gridpoints;

           x = reshape(vol.interpolant(grxy),diff(me.axdim([1 3; 2 4]))+1)';
%            trfunm = @(x)x';
%                    me.rtrmat(me);
            irg =vv.intensity_range;
             x = (x-irg(1))/diff(irg);
             x(x>1) = 1;
             x(x<0) = 0;

            me.ploth = me.ploth(ishandle(me.ploth));
%             plh = ismember(vol.ploth,get(me.h,'children'));
              
%             set(me.ploth,'cdata',x(:,:,[1 1 1]),'visible','on','xdata',[0 size(x,2)-1],'ydata',[0 size(x,1)-1]);                          
            set(me.ploth,'cdata',x(:,:,[1 1 1]),'visible','on','xdata',me.axdim([1 2]),'ydata',me.axdim([3 4]));                          
             
            %%% Now plot cross hairs whos arms depict the intersection of
            %%% all other view planes.
              
             crk = me.handles.crosshairs;   
             crk = crk(ishandle(crk));
             if me.showcrossHere && get(vv.showcrossh,'value')
                 
                 other_axes = vv.plotax(vv.plotax~=me & isopen(vv.plotax));
                 
                 crp = vv.current_point;
                 nv = me.normvec;
                 for kk = 1:length(other_axes)
                    oax = other_axes(kk);
                    if oax.showcrossThere
%                         axl=axis(oax.h);
                        axl = [get(oax.h,'XLim'),get(oax.h,'YLim')];
                        cpax = oax.vol2ax(vv.current_point)';
                        xc = [axl(1:2)',cpax([ 2 2]);cpax([1 1]),axl(3:4)'];
                        nv2 = oax.normvec(1:3);
                        vx = oax.ax2vol(xc);
                        q=cross(nv,nv2);
                        q=q./norm(q);
                        qcx= ((vx-repmat(crp,size(vx,1),1))*q');
                        [~,mni] = min(qcx);
                        [~,mxi] = max(qcx);
                        vxpr = qcx([mni mxi],:)*q+repmat(crp,2,1); 
                        crlim = me.vol2ax(vxpr);
                        set(crk(kk),'xdata',crlim(:,1),'ydata',crlim(:,2),'visible','on','zdata',[.5 .5],'color',oax.crosscol,'linestyle',oax.crossprop)
                 
                    else
                         set(crk(kk),'visible','off')
                    end
                 end
             else
                 set(crk(:),'visible','off')
             end
             
             if isprop(me,'scalebar') && me.scalebar.length>0
                
                 dd = @(x)sqrt(sum(x.^2,2));
                 plotlen = dd([me.scalebar.length 0 0 ]*me.parent.transforms(1).trmat(1:3,1:3)^-1*me.trmat(1:3,1:3));
%                  axdim = axis(me.h);
                 axdim = [get(me.h,'XLim'),get(me.h,'YLim')];
                 if ishandle(me.scalebar.handle)
                     delete(me.scalebar.handle)
                 end
                 plotx = diff(axdim(1:2))*me.scalebar.position(1)+axdim(1);
                 ploty = diff(axdim(3:4))*me.scalebar.position(2)+axdim(3);
                 
                 me.scalebar.handle = rectangle('position',[plotx,ploty,me.scalebar.length.*[ 1 me.scalebar.width]],...
                                                'facecolor',me.scalebar.color,'edgecolor','none','parent',me.h);
                 
             end
            if get(vv.fixSisterAx,'value') && updatesis
                 sisax = me.sisobj;
                 if ~isempty(sisax)
                     sisax(~isvalid(sisax))=[];
                     sisax(~ishandle([sisax.h])) = [];
                     for kk = 1:length(sisax)
                            sisax(kk).plotupdate(false);
                            axis(sisax(kk).h,axis(me.h))
                     end
                     me.sisobj=sisax;
                 end
            end
            
    end
    %%%
    function a = get.sisters(me)
       a = me.sisobj; 
    end
    
    function  set.sisters(me,a)
       
       if isempty(a)
           me.sisobj = me.sisobj([]);
       else
           me.makesis(a(isvalid(a)));
       end   
       
    end
      
   
    function makesis(me,a,reciprocal)
      if nargin < 3
            reciprocal =true;
      end
       a = a(a~=me);
       a = a(isa(a,'plotax')& isvalid(a));
       me.sisobj = unique([me.sisobj(:)',a]);    
       if reciprocal
           for k = 1:length(a)
              a(k).makesis(me,false);        
           end  
       end
    end
    
    %%%
    
    function reset(me)

        vol = me.parent.current.volumes(1); 
        a = @(x)[[x;x],kron([0 1]',ones(size(x,1),1))];         

        dperm = [1 2 3];
        sz = size(vol.image.Data);
        corners = a(a([0 1]'))*diag(sz(dperm)-1);
        %          dci = chooseperm(size(corners,1),2);
        corners(:,4)=1;

        u = me.handles.control;
        nv = [me.normvec 0];

        zlim = minmax(corners*nv');

        set(u,'Min',zlim(1),'Max',zlim(2),'SliderStep',[1/diff(zlim) .1],'value',me.parent.current_point*nv(1:3)'); 
        me.plotupdate;
        axis(me.h,'image','xy');
        if isempty(me.axdim) || max(diff(me.axdim))<2
%              me.axdim=axis(me.h);
              axdim = [get(me.h,'XLim'),get(me.h,'YLim')];
        end

    end
    
    
   end  %%% End of Methods block
 
end