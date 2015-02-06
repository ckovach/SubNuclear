
classdef volumeview < handle
    
%
%  The volumeview class allows viewing slices through image volumes with 
%  mesh and point objects superimposed. Through the associated graphical
%  interface it also allows for the selection and editing of mesh and point
%  data.
%
%  v = volumeview initiates a new class and prompts the user to give a
%  volume file.
%
%  v.addmesh(Tr) adds the mesh contained in a TriRep or mesh object.
%
%  v.sisters = v2 ties axis properties between two volumeview objects to
%               each other.
%         
%  v2 = volumeview(v) creates a new volumeview object with the same attributes
%                   as v (i.e. a copy of v) and initializes the graphical
%                   interface. To bring up a GUI for a volumeview object which
%                   has been loaded from a file, use v = volumeview(v).
% 
%  Note that because volumeview is a handle class object 
%         v2 = v
%  does not create a copy of v and v2 is only another reference to the same
%  object.
% 
% PROPERTIES
%
% current_point :  Current point in voxel coordinates
% current_coord :  Current point in transformed coordinates
% sisters : Sister volumeview objects.
% intensity_range : Scaling of the image plot.
% current_object : Currently active object(s) (point, mesh or transform)
% meshes : Array of mesh objects representing 3D surfaces in the image space. 
% points : Array of point objects representing points in the image space.
% volumes : Volume (currently non-array)/
% transforms : Array of transform objects representing affine
%               transformations of voxel coordinates.
% info : General information about the volumeview object.
% path :  Path to file.
% fig :  Figure handle for main GUI.
% plotax : axis object.
% current : structure with whos fields give the most recently selected 
%            objects of each category as current.(object_type), where
%            (object_type) is points, meshes, transforms, or volumes.
%
% METHODS
% 
% addtransform : Add a transform.                     
% addpoint : Add a point.
% addmesh : Add a mesh.
% resetaxis : Reset axis plots.           
% makeslice :  Returns a slice of arbitrary orientation.
%
% plotupdate : update plots.       
% makesis : make an object a sister.
%
% See also MESHES POINTS TRANSFORMS

% C Kovach 2013

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

 properties

     meshes ;
     points ;
     volumes;
     transforms ;
     parent ;

 end
 properties (Dependent = true)
     current_point;
%      Vol = [];
         sisters ;
         intensity_range;
     current_object;
     current_coord ;
     
    objectiter; 
 end
 properties (SetAccess = private)
%      meshes =meshes;
%      points =points;
%      volumes =volumes;
%      transforms =transforms;
%      info =info;
       modules;
      modobjects;
     plotax;
     info ;
     path = '';
     activeax = [];
     activemesh = [];
     activeimage = [];
     fig = [];
%      axes = [];
     sisobj = [];
     current = [];
     isplotting = false; %True while plotting has not finished.
 end
  properties (SetAccess = private, Hidden = true)
%     start_transposed = false; %images initially transposed
    currp =[0 0 0];
    sliders ;
     objiter = 0;
    pointbox ;
    coordbox ;
%      imhandles ;
     bdownpt=[];
     buppt=[];
       meshadjust= [];
    moving_mesh = false;
    meshlist ;
    ptlist ;
    vollist ;
    transflist ;
    modulelist ;
    annoh ;
    addpt ;
    showcrossh ;
    unitlabel ;
%     objectiter = 0;
%     current_volume ;
%     crosshairs ;
    croj = struct('field','info','id',0);
    annohtitle=[];
    legendh ;
    showlegendh ;
    fixSisterAx ;
    extrafigs ;
    object_types = {'meshes','points','volumes','info','transforms','plotax'};%,'lines'};
    currtypes= [];
    currlist ;
    blowh ;
    suckh ;
    plprop ;%
% prop = [];
  end
  
 methods
     
     function     me = volumeview(varargin)
 
         if nargin == 1 && isempty(varargin{1}) 
            me = me([]);
            return
         end
         axsz = [.35 .4];
         axpos = [.05 .55 axsz
                 .55 .55 axsz
                 .3 .05 axsz];
        
%          me.image = medimage(varargin{:}); % load volume
         objects = me.object_types;
         for k = 1:length(objects)
              me.(objects{k}) = feval(str2func(objects{k}));
             switch objects{k}
                 case {'info'}
            otherwise                    
                     me.(objects{k})= me.(objects{k}).empty;
             end
             me.currtypes.(objects{k})= [];

         end        
        me.fig = figure('NumberTitle','off');        
%          set(me.fig,'WindowButtonMotionFcn',@(a,b)me.wmvfn(a,b),...
%            'WindowButtonUpFcn',@(a,b)me.figButtonUp(a,b),'ButtonDownFcn',@(a,b)me.figButtonDownFcn)
        me.addfig(me.fig);
      
          if nargin > 0 && isa(varargin{1},'volumeview')
              vw = varargin{1};
              props= setdiff(properties(vw),{'fig','plotax','sisobj','sisters'});
              for i =1:length(props)
                  try 
                      me.(props{i}) = vw.(props{i});                  
                  end
              end
              for i = 1:length(objects)
                  for k = 1:length(me.(objects{i}))
                     me.(objects{i})(k).objectid= me.objectiter;
                  end
                  me.info.objectid = 0;
              end

              lbl = me.info.label;
              me.plotax(~ishandle([me.plotax.h])) = [];
     
              for i = 1:length(me.object_types)
                  objs = me.(me.object_types{i});
                  if isa(objs,'plotobj')
                      for k = 1:length(objs)
                          objs(k).ploth = objs(k).ploth(ishandle(objs(k).ploth));
                      end
                  end
              end
              me.sisobj = [];
          elseif nargin >0 && ischar(varargin{1})
             lbl = sprintf('SubNuclear: %s',varargin{1});
             me.info.label = lbl;
         else
             lbl = 'SubNuclear';
          end
         set(me.fig,'name',lbl)
          %%%% Set up GUI   
         %%% Current point in units and coordinates
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.3 .4], .2,.03],'string','Current Point','fontsize',12);
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.25 .33], .02,.04],'string','X');
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.25 .29], .02,.04],'string','Y');
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.25 .25], .02,.04],'string','Z');
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.24 .38], .03,.02],'string','vox');
        me.unitlabel = uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.21 .38], .04,.02],'string','units');
        pt(1) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.23 .33], .03,.04],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.ptboxupdate(a,b));
        pt(2) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.23 .29], .03,.04],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.ptboxupdate(a,b));
        pt(3) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.23 .25], .03,.04],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.ptboxupdate(a,b));
        cpt(1) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.2 .33], .03,.04],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.coordboxupdate(a,b));
        cpt(2) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.2 .29], .03,.04],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.coordboxupdate(a,b));
        cpt(3) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.2 .25], .03,.04],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.coordboxupdate(a,b));
        
        %%% Radio button to show cross hairs and add points
        me.showcrossh = uicontrol('style','radiobutton','units','normalized','position',[axpos(3,1:2)+[-.17 .32],.05,.04],'string','+','Callback',@(a,b)me.plotupdate(),'fontsize',16);
        me.addpt = uicontrol('style','pushbutton','units','normalized','position',[axpos(3,1:2)+[-.17 .27], .05,.04],'string','Add','Callback',@(a,b)me.addpoint);        
        me.pointbox = pt;   
        me.coordbox = cpt;   
        %Push button to contract and expand meshes
        me.suckh = uicontrol('style','pushbutton','units','normalized','position',[axpos(3,1:2)+[-.10 .32], .05,.04],...
                              'string','Contract','Callback',@(a,b)me.mesheditCback(a));         
        me.blowh = uicontrol('style','pushbutton','units','normalized','position',[axpos(3,1:2)+[-.10 .27], .05,.04],...
                              'string','Expand','value',0,'Callback',@(a,b)me.mesheditCback(a));
        % Context menu for list boxes
        ucm = uicontextmenu;
        uimenu(ucm,'label','Goto','callback',@(src,evnt)me.gotoobj(src,evnt));
        uimenu(ucm,'label','Toggle hide','callback',@(src,evnt)me.hideobj(src,evnt));
        uimenu(ucm,'label','Label','callback',@(src,evnt)me.renameobj(src,evnt));
        uimenu(ucm,'label','Load','callback',@(src,evnt)me.loadobject(src,evnt));
        uimenu(ucm,'label','Add','callback',@(src,evnt)me.addobject(src,evnt));
        uimenu(ucm,'label','Copy','callback',@(src,evnt)me.copyobject(src,evnt));
        uimenu(ucm,'label','Write','callback',@(src,evnt)me.writemesh());
        uimenu(ucm,'label','Help','callback',@(src,evnt)me.help_obj());
        uimenu(ucm,'label','-------','callback',@(src,evnt)[]);
        uimenu(ucm,'label','Delete','callback',@(src,evnt)me.delobj(src,evnt));
        
        ucmodhelp = uicontextmenu;
        uimenu(ucmodhelp,'label','Help','callback',@(src,evnt)me.modhelp(src,evnt));
      
        %%% radio button to enable mesh editing
        adj = uicontrol('style','radiobutton','units','normalized','position',[axpos(3,1:2)+[.4 .35], .1,.03],'string','Adjust Mesh');
        me.meshadjust = adj;
        me.fixSisterAx = uicontrol('style','radiobutton','units','normalized','position',[axpos(3,1:2)+[.4 .32], .1,.03],...
                                  'string','Link axes to sisters','value',1,'Callback',@(a,b)me.plotupdate());
        
        %%% Edit boxes to control plotting properties of selected object
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[.51 .35], .03,.03],'string','color');
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[.54 .35], .03,.03],'string','style');
        uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[.57 .35], .05,.03],'string','plot args');
%         uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[.65 .35], .05,.03],'string','Marker');
        me.plprop(1) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[.51 .32], .03,.03]);
        me.plprop(2) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[.54 .32], .03,.03]);
        me.plprop(3) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[.57 .32], .05,.03]);
        me.showlegendh = uicontrol('style','checkbox','units','normalized','position',[axpos(3,1:2)+[.62 .32], .05,.03],'string','legend','value',0,'callback',@(a,b)me.legendcheck());
%         me.plprop(3) = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[.65 .32], .05,.03]);
        set(me.plprop,'Callback',@(a,b,c) me.updatePlotProps(a));
        me.meshlist  = uicontrol('style','list','units','normalized','position',[axpos(3,1:2)+[.39 0], .1,.3],'string','Meshes',...
                                 'Callback',@(a,b)me.listcallback(a,b),'backgroundcolor',[1 1 1],'uicontextmenu',ucm,'max',2);
        
        me.ptlist  = uicontrol('style','list','units','normalized','position',[axpos(3,1:2)+[.49 0], .1,.3],'string','Points',...
                                'Callback',@(a,b)me.listcallback(a,b),'backgroundcolor',[1 1 1],'uicontextmenu',ucm,'max',2);
%         me.vollist  = uicontrol('style','list','units','normalized','position',[axpos(3,1:2)+[.59 0], .1,.3],'string','Volumes',...
%                                  'Callback',@(a,b)me.listcallback(a,b),'backgroundcolor',[1 1 1],'uicontextmenu',ucm,'max',2);
        me.transflist = uicontrol('style','list','units','normalized','position',[axpos(3,1:2)+[.59 0], .1,.15],'string','Transforms',...
                                 'Callback',@(a,b)me.listcallback(a,b),'backgroundcolor',[1 1 1],'uicontextmenu',ucm,'max',2);    
        me.modulelist = uicontrol('style','list','units','normalized','position',[axpos(3,1:2)+[.59 .15], .1,.15],'string','MODULES',...
                                 'Callback',@(a,b)me.listcallback(a,b),'backgroundcolor',[1 1 1],'uicontextmenu',ucmodhelp,'max',2);    
    
        %%% Annotation listbox
        me.annoh = uicontrol('style','edit','units','normalized','position',[axpos(3,1:2)+[-.25 0], .25,.2],'string','','max',2,'horizontalalignment','left',...
                              'Callback',@(a,b) me.AnnoCallback(a,b),'backgroundcolor',[1 .949 .86]);
        me.annohtitle = uicontrol('style','text','units','normalized','position',[axpos(3,1:2)+[-.25 .2], .25,.03],'string',['Annotation: ',me.info.label],'fontsize',12);
        %%%%%%%
        
        sm = what('SubModules');
        mods = dir(fullfile(sm(1).path,'*.m'));
        me.modules = setdiff(arrayfun(@(x)x.name(1:end-2),mods,'uniformoutput',0),'module');
          mdo = me.modules(:)';
        mdo(2,:) = {[]};
        me.modobjects = struct(mdo{:});
                 % Begin by adding a volume
         if nargin < 1 || ~isa(varargin{1},'volumeview')
             me.addvolume(varargin{:});
         end
        for k = 1:3
            
             me.addaxis(axes('position',axpos(k,:),'units','normalized','parent',me.fig)); %#ok<LAXES>
%              sz = size(me.Vol);
%              nv = double(1:3==k);
            e = eye(3);
%             T = [e(setdiff(1:3,k),:),[1 1]']';
            T = e(setdiff(1:3,k),:)';
            T(:,3) = sum(T,2) ==0;
            T(4,4) = 1;
%             me.plotax(k).axdim = sz*T(:,1:2);
            stdtr = me.transforms(1)*transforms('trmat',diag([-1 1 1 1])*T); %#ok<*CPROP>
%             stdtr = transforms('trmat',diag([-1 1 1 1])*T); %#ok<*CPROP>
            me.plotax(k).Transform = stdtr;
%             me.plotax(k).normvec = stdtr.tr(nv);
% %             me.plotax(k).axvec= e(setdiff(1:3,k),:);
%             me.plotax(k).trmat= T;
            if ~isempty(me.Vol)
                me.plotax(k).setvolmat(size(me.Vol)/2);
            end   
        end
        if isempty(me.transforms)
            me.addtransform(eye(4),'Identity');
        end
        me.current.transforms = me.transforms(1);
%         me.sliders = u;
       
      
         me.resetaxis();
%         me.plotupdate;
          
     end
     %%%

     function a=copy(me) % Create a copy of the object 
         a=volumeview(me);
         a.resetaxis;
     end

     function a = Vol(me,varargin)         
         a = me.current.volumes(1).image.Data(varargin{:});
     end
     
     %%%%
     
     function set.intensity_range(me,a)
         
         cv = me.current.volumes;
         for i = 1:length(cv)
             cv(i).intensity_range = a;
         end
         me.plotupdate();
         
     end
     %%%
     
     function set.current_object(me,obj)
         
         if isempty(obj), return, end
         hh = cellfun(@(fld)[me.(fld).objectid],me.object_types,'uniformoutput',0);
         fld = me.object_types{cellfun(@(x)any(ismember([obj.objectid],x)),hh)};
         
         obji = ismember([me.(fld).objectid],[obj.objectid]);
         
%          me.croj.field = fld;
%          me.croj.objectid = [me.(fld)(obji).objectid];
%          me.current.(fld)=[me.(fld)(obji).objectid];
         fld = class(obj);
         me.current.object= obj;
         me.current.(fld)= obj;
         
         switch fld
             case 'volumes'
                 lst = me.vollist;
                 case {'points','lines'}
                 lst = me.ptlist;
                 
             case 'meshes'
                 lst = me.meshlist;
             case 'transforms'
                 lst = me.transflist;
             otherwise 
                 return
         end
         set(lst,'value',find(obji)+1);
         
         if isa(obj,'plotobj')
             plc = obj(1).plotcolor;
             if isnumeric(plc)
                 plc = strcat('[',num2str(plc),']');
             end
            set(me.plprop(1),'string',plc);
            set(me.plprop(2),'string',obj(1).linestyle);
            plarg = obj(1).plotargs;
            isn = cellfun(@isnumeric,plarg);
            plarg(isn) = cellfun(@(x)mat2str(x),plarg(isn),'uniformoutput',false);
            str = sprintf('%s,',plarg{:});
            set(me.plprop(3),'string',str(1:end-1));
            
%             set(me.plprop(3),'string',obj.marker);
            
         end
         
     end
     
     function a = get.current_object(me)
         
%          fld = me.croj.field;         
%          obji = ismember([me.croj.objectid],me.(fld).objectid);
         
%          a = me.(fld)(obji);
            a = me.current.object;
     end
     function a = get.intensity_range(me)
         cv = me.current.volumes;
%          imi = ismember([me.volumes],cv);
         a = cv.intensity_range;
     end
     %%%%
     function a = get.current_point(me)         
         a = me.currp;
     end
     %%%
     function set.current_point(me,a)
        me.currp = a;
        for i= find(isopen(me.plotax))
%             ctrl = me.plotax(i).handles.control;
             me.plotax(i).setvolmat(a);

%             set(ctrl,'value',a*nv')
            
        end
        try
            me.plotupdate;
        catch err
            warning('Caught error while updating plot:\n\t%s\n',err.message)
            me.cleanup;
            me.plotupdate;
        end
        if get(me.fixSisterAx,'value')
         for k = 1:length(me.sisters)
             sis = me.sisters(k);
             sis.currp= me.current_point;
             sis.plotupdate(sis.plotax,false);
          end
        end
     end
     
     %%%
     
     function a = get.objectiter(me)
        % To ensure a unique id for each object, the object iterator
        % increments everytime it is called.
         me.objiter = me.objiter+1;
         a = me.objiter;
         
     end
     function set.objectiter(me,a)
        
        error('Object iterator cannot be reset.')
     end
     %%%
     
     %%%
%      function ptvol = ax2vol(me,ax,ptax)
%                
%         ptvol = ax.ax2vol(ptax);
%          
%      end

     %%%
     
     function addtransform(me,T,label,varargin)
         
         
          ntrf = length(me.transforms);
          if nargin < 2
              [fn,pth] = uigetfile;
              T = fullfile(pth,fn);
          end
          if nargin <3 && ~isa(T,'transforms')
             label = sprintf('Transform %i',ntrf+1);
             label = inputdlg('Enter Label','Label transform',1,{label});           
          elseif isa(T,'transforms')
              label = {T.label};
          elseif ~iscell(label)
              label = {label};
          end
        if ischar(T)
            T = textread(T); %#ok<REMFF1>
            if size(T,2)==5
                T = T';
            end
            T = T(1:4,1:4);
        end
        if ~isa(T,'transforms')
            me.transforms(ntrf+1) =transforms('label',label{1},'trmat',T,'file','','objectid',me.objectiter,'notes',''); %#ok<CPROP>
        else
            T.objectid = me.objectiter;
            T.label = label{1};
            me.transforms(end+1) = T;
        end
            me.current_object = me.transforms(ntrf+1);
            me.updatelists();
%         set(me.transflist,'value',ntrf+1)
     end
     %%%
     function loadobject(me,src,evnt) %#ok<*INUSD>
          switch me.currlist
              
              case me.meshlist
                  me.addmesh();
              case me.transflist
%                   me.loadtransform(); %Not implemented yet
              case me.ptlist
                  [fn,pth]=uigetfile('*.*');
                  me.addpoint('',fullfile(pth,fn));
          end
     
     end
     function addobject(me,src,evnt)
          switch me.currlist
              case me.meshlist
                  me.addmesh();
              case me.ptlist
                  me.addpoint();
              case me.transflist
                  maketransform(me,me);          
          end
     end
     function copyobject(me,src,evnt)
        co = me.current_object;
        fld = class(co);
        switch fld
            case 'meshes'
                
                me.addmesh(co);
                newm=me.meshes(end);
                newm.label = ['Copy of ',newm.label];
            case 'points'
                for i = 1:length(co)
                    me.addpoint(co(i));
                    me.points(end).label = ['Copy of ',co(i).label];
                end
            case 'volumes'
        end
     end
     %%
     function addmesh(me,trirep,label,varargin)
         nmesh = length(me.meshes);
          if nargin <3 
            label = {sprintf('mesh %i',nmesh+1)};
          elseif ~iscell(label)
              label = {label};
          end
         if nargin < 2 ||  ~isa(trirep,'TriRep') && ~isa(trirep,'meshes')

%              file = '';
             if nargin < 2 || isempty(trirep)
                [fn,pth] = uigetfile({'*.vtk;*.bmv','Mesh Files';'*.*','All Files'}, 'Load Mesh','multiselect','on');
             elseif ~isa(trirep,'mshes')
                 [pth,fn,ext]= fileparts(trirep);
                 fn = [fn,ext];
             end
             if ~iscell(fn)
                 fn = {fn};
             end
             for kk = 1:length(fn)
                 nmesh = length(me.meshes);
                 if nargin <3 
                    label = {sprintf('mesh %i',nmesh+1)};
                 else
                     label = label(kk);
                 end
                file = fullfile(pth,fn{kk});
                vtk = readvtk(file);
                trirep = TriRep(vtk.tri,vtk.vert);
    %              me.objectiter = me.objectiter+1;
                 me.meshes(nmesh+1).label = label{1};
                 me.meshes(nmesh+1).trirep = trirep;
                 me.meshes(nmesh+1).file= file;
                  me.meshes(nmesh+1).show = true;
                  me.meshes(nmesh+1).objectid = me.objectiter;
                  me.meshes(nmesh+1).notes = '';
                 me.current_object = me.meshes(nmesh+1);
                for k = 1:length(me.plotax)
                    plh = plot(me.plotax(k).h,0,0,'.',varargin{:});
                    me.meshes(nmesh+1).ploth(k) =plh;
                    set(plh,'ButtonDownFcn',@(a,b)me.meshButtonDown(a,b))
                end
                if nargin < 3
                    me.renameobj(fn{kk});
                end
             end
         elseif isa(trirep,'meshes')
            for i = 1:length(trirep)
                props = properties(trirep);
                newm = meshes;
                for k = 1:length(props)
                    try %#ok<*TRYNC>
                        newm.(props{k}) = trirep(i).(props{k});
                    end
                end
%                 me.objectiter = me.objectiter+1;
                  
                  newm(end).objectid = me.objectiter;
            
                  me.meshes(end+1) = newm; 
               for k = 1:length(me.plotax)
                    plh = plot(me.plotax(k).h,0,0,'.',varargin{:});
                    me.meshes(end).ploth(k) =plh;
                    set(plh,'ButtonDownFcn',@(a,b)me.meshButtonDown(a,b))
               end
         
            end
         else
                 me.meshes(end+1).trirep = trirep;
%                 me.objectiter = me.objectiter+1;
                me.meshes(end).objectid = me.objectiter;
                me.meshes(end).label = label{1};
                me.meshes(end).show = true;
               for k = 1:length(me.plotax)
                    plh = plot(me.plotax(k).h,0,0,'-',varargin{:});
                    me.meshes(end).ploth(k) =plh;
                    set(plh,'ButtonDownFcn',@(a,b)me.meshButtonDown(a,b))
               end

%              me.meshes(end)
%              error('Input must be a trirep object, a file name, or a mesh structure');
         end
         me.current_object = me.meshes(end);
        me.annotationUpdate();
        me.plotupdate();
        me.updatelists();
             
     end
     %%%

     %%%
     function addvolume(me,varargin)
        
%          if nargin < 3 
%             vol =[];
%          elseif nargin >= 3
%              vol = varargin{1};
%          end
         if nargin < 2 || ~ischar(varargin{1})
             label = '';
         else
             label = varargin{1};
       %      varargin(1) = [];
         end
%          if isempty(vol)
%               vol = medimage(varargin{:});
%          end
%          if ~isa(vol,'medimage')
            vol = medimage(varargin{:}); 
%          end
%          vol.reorient2std;
         nvol = length(me.volumes);
        % voli = nvol+1;
         
         sz4 = size(vol.Data,4);
         if sz4 >1
             lblstr = '%s %i';
             for k = 1:sz4;
                 vv(k) = medimage(vol);
                 vv(k).Data = vol.Data(:,:,:,k);
             end
             vol = vv;
         else
             lblstr = '%s';
         end
        if nargin <2 || isempty(label)
            label = 'Volume';
            label = inputdlg('Enter Label','Label volume',1,{label});           
            label = label{1};
        end     
%          vol.reorient2std;
         for k = 1:sz4
              voli = nvol+k;
            
            me.volumes(voli).label = sprintf(lblstr,label,k);
         
             me.volumes(voli).image= vol(k);
             me.volumes(voli).file= '';
              me.volumes(voli).show = true;
    %           me.objectiter=me.objectiter+1;
              me.volumes(voli).objectid = me.objectiter;
              me.volumes(voli).notes = '';
            me.volumes(voli).intensity_range = minmax(vol(k).Data(:));
            me.volumes(voli).parent = me;

            me.current_object =  me.volumes(nvol+1);
    %         me.current.volumes =  me.volumes(voli);

            for k = 1:length(me.plotax)

                            me.plotax(k).ploth(voli) = image(nan,'parent',me.plotax(k).h); %#ok<CPROP>

                    set( me.plotax(k).ploth(voli),'ButtonDownFcn',@(a,b)me.plotax(k).axisupdate(a,b))
                    me.plotax(k).reset;
           end
    %         me.imhandles(voli,:) =  me.volumes(voli).ploth;

            tr = [me.volumes(voli).image.vox2mm; 0 0 0 1];
            if ~isempty(tr)
                me.addtransform(tr',sprintf('vol%i mm',voli));
            end
         end
          me.updatelists();
          tr = me.volumes(voli).tr2std;
          tr.label = sprintf('vol%i std',voli);
          me.addtransform(tr);
          me.annotationUpdate();
         me.plotupdate();
     end
     %%%
     function cleanup(me,cleansis)
         if nargin  < 2
             cleansis = true;
         end
         me.plotax(~isvalid(me.plotax))=[];
         me.plotax(~ishandle([me.plotax.h]))=[];
         for i = 1:length(me.object_types)
             fld = me.object_types{i};
             if isa(me.(fld),'plotobj') 
                 me.(fld)(~isvalid(me.(fld))) = [];
                 for k = 1:length(me.(fld))                    
                     me.(fld)(k).ploth(me.(fld)(k).ploth==0 | ~ishandle(me.(fld)(k).ploth)) = [];
                     me.(fld)(k).ploth(length(me.plotax)+1:end)=[];
                 end
             end
         end
%          me.plotax(~ismember([me.plotax.h],[me.volumes.ploth]))=[];
         if ~isempty(me.sisobj)
            me.sisobj(~isvalid(me.sisobj)) = [];
            if cleansis
                for i = 1:length(me.sisobj)
                   me.sisobj(i).cleanup(false); 
                end
            end
%             arrayfun(@(x)x.cleanup(false),me.sisobj);
         end
     end
     
     %%%%%%%%
     function resetaxis(me,ax)
         %%% Update axis with new trdat and return to max scale
         if nargin < 2
             ax = me.plotax;
         end
         if isnumeric(ax)
             ax = me.plotax(ax);
         end
         me.cleanup();
         
      
%            dc =corners(dci(:,2),:)-corners(dci(:,1),:);
%           plh = cat(1,me.volumes.ploth);
         for k = 1:length(ax)
           set(ax(k).ploth,'ButtonDownFcn',@(src,evnt)ax(k).axisupdate(src,evnt))
              ax(k).reset;
         end
        axis([ax.h],'equal');
         lsts = [ me.meshlist ,me.ptlist,me.vollist ,me.transflist ,me.modulelist] ;
         set(lsts,'value',1)
%          set([me.volumes.ploth],'ButtonDownFcn',@(src,evnt)me.axisupdate(src,evnt))
     end
     
     %%%%%%%
     function addpoint(me,label,varargin)
          
         trfun = @(x)x;
         plotargs = {};
         if nargin >1 && isa(label,'points')            
             pt = label;
             label = pt.label;
         elseif nargin < 3 || isempty(varargin{1})
             pt = me.current_point;
         else
             pt = varargin{1};
             varargin(1) = [];
         end
         
         i = 1;
         while i <= length(varargin)
             switch varargin{i}
                 case 'transform'
                     if isnumeric(varargin{i+1})
                         transform = varargin{i+1};
                         transform(4,1:4) = [0 0 0 1];
                         trfun = @(x)[x 1]*transform;
                     elseif isa(varargin{i+1},'function_handle');
                         trfun = varargin{i+1};
                     elseif isa(varargin{i+1},'transforms')
                         trfun = varargin{i+1}.tr;
                     end
                     i = i+1;
                 case 'plotargs'
                     plotargs = varargin{i+1};
                     i = i+1;
                     
                 otherwise 
                     error('Unrecognized option %s',varargin{i})
             end
             i=i+1;
         end
         
         npt = length(me.points);
         if nargin <2 || isempty(label)
            label = sprintf('Point %i',npt+1);
            label = inputdlg('Enter Label','Label point',1,{label});           
            label = label{1};
         end      
         
         
         if ~isa(pt,'points')
             if size(pt,1)>1
                 label = [label,' %i'];
             end
             if ischar(pt)
                 pt = dlmread(pt);
             end
             for k = 1:size(pt,1)
                 me.points(npt+k).label = sprintf(label,k);
                 me.points(npt+k).coord= trfun(pt(k,:));
                 me.points(npt+k).file= '';
                  me.points(npt+k).show = true;
        %           me.objectiter=me.objectiter+1;
                  me.points(npt+k).objectid = me.objectiter;
                  me.points(npt+k).notes = '';
             end
            
              newpt = me.points(npt+1:end);
         else
            
             getprops = setdiff(properties(pt),{'objectid','coord'});
             for i = 1:length(pt)
                 me.points(npt+i).objectid = me.objectiter;
                  me.points(npt+i).coord = trfun(pt(i).coord);
                 for k = 1:length(getprops)
                    me.points(npt+i).(getprops{k}) = pt(i).(getprops{k});
                    
                 end
%                  if nargin < 2
%                      me.points(npt+i).label = ['Copy of ',me.points(npt+i).label ];
%                  end
             end
             newpt =  me.points(npt+1:end);
            
         end
         me.current_object = me.points(npt+1);
         for i = 1:length(newpt)
            for k = 1:length(me.plotax)            
                plh = plot(me.plotax(k).h,0,0,'.',plotargs{:});
                newpt(i).ploth(k) =plh;
                set(plh,'ButtonDownFcn',@(a,b)me.meshButtonDown(a,b))
            end
         end
        me.annotationUpdate();
        me.plotupdate();
           me.updatelists();
     end
    
     %%%
     function annotationUpdate(me,src,evnt)
        cro =  me.current_object(1);
        fld = class(cro);
        indx = ismember([me.(fld).objectid] , [cro.objectid]);
        
        obj = me.(fld)(indx);
        if isempty(obj),return, end
        set(me.annohtitle,'string',sprintf('Annotation: %s',obj(1).label));
        set(me.annoh,'string',obj(1).notes);
        
     end
     %%%
     function figButtonDownFcn(me,src,evnt)
         
        me.current_object=me.info;
        
     end
         
     %%% Make another volumeview object a sister
     function makesis(me,sis,reciprocal)
         if nargin < 3 
             reciprocal = true;
         end
         
        if ~isempty(sis)
         for i = 1:length(sis)
             nsisobj(i) = sis(i); %#ok<*AGROW>
             sis(i).current_point = me.current_point;
             for k = 1:3
                 me.plotax(k).sisters = [me.plotax(k).sisters,sis(i).plotax(k)];
             end
%              sis(i).plotupdate();
             if reciprocal
                 sis(i).makesis(me,false);
             end
         end
        else
            nsisobj = [];
        end
         me.sisobj = [me.sisobj,nsisobj];
         me.plotupdate();
     end        
     %%%
     function set.sisters(me,sis)
         if isempty(sis)
             me.sisobj = me.sisobj([]);
         else
             me.makesis(sis(sis~=me));
         end
     end
     %%%
     function b=get.sisters(me)
         b=me.sisobj;
     end
     

    %%%
    function addfig(me,nfig)
        
        me.extrafigs(end+1) = nfig;
      set(nfig,'DeleteFcn',@(a,b)me.rmfig(a,b),'WindowButtonMotionFcn',@(a,b)me.wmvfn(a,b),...
       'WindowButtonUpFcn',@(a,b)me.figButtonUp(a,b),'ButtonDownFcn',@(a,b)me.figButtonDownFcn);

    end
    %%%
    
    function rmfig(me,src,evnt)
        me.fig(me.fig == src) = [];
        me.extrafigs(me.extrafigs == src) = [];
    end
    %%%
    function axout = addaxis(me,axh,addsis)
        
        % Add a new axis to the volumeview object
        %
        % Use: vv.addaxis(axis_handle)
        
        if nargin < 3
            addsis = true;
        end
        
        ax = plotax(me,axh);
%         ax.parent = me;
        
        me.plotax(end+1) = ax;

        axi = length(me.plotax);
        hold(ax.h,'on');
        flds = me.object_types;
        for i = 1:length(flds)
            fld =flds{i};
            switch fld
                case 'meshes'
                    mk = '.';
                case 'points'
                    mk = '+';
            end
                    
            for k = 1:length(me.(fld))
                switch fld
                    case {'meshes','points'}
                        me.(fld)(k).ploth(end+1) = plot(axh,nan(1),nan(1),mk);
                        set(me.(fld)(k).ploth(end),'ButtonDownFcn',@(a,b)me.meshButtonDown(a,b))
                    case {'volumes'}
                        im = image(nan,'parent',axh); %#ok<CPROP>
                        me.plotax(axi).ploth(k)=im;
                        set(im,'visible','off','ButtonDownFcn',@(a,b)ax.axisupdate(a,b))
                    case {'info','transforms','plotax'}
                    otherwise
                        error('Unrecognized field. Not sure what to do...')

                end
                        
            end
        end
      
        axout = ax;
        if addsis
            for i = 1:length(me.sisters)
               fg = figure;
               axs = me.sisters(i).addaxis(axes('parent',fg),false); %#ok<LAXES>
               axs.basetrmat = ax.basetrmat;
               axs.Transform = ax.Transform;
%                axs.normvec = ax.normvec;
               axs.axdim = ax.axdim;
               axis(axs.h,'equal');
               axis(axs.h,axis(ax.h));
               ax.sisters = [ax.sisters,axs];
               axs.sisters = [axs.sisters,ax];
%                 axs.parent.plotupdate();
                axout = [axout,axs.sisters];
            end
        end
    end
    
   
   
      %%%%%
    function a = get.current_coord(me)
        crT = me.current.transforms;
        a = crT(1).tr(me.current_point);        
    end
     %%%
    function set.current_coord(me,a)
      
        me.current_point = crT.itr(a);
        
    end
    function mesheditCback(me,src,evnt)
        
        switch src
            case me.suckh
                val  = 1;

            case me.blowh
                val = -1;

        end
        me.meshedit(val);
    end
    
    
 end
    methods (Hidden = true)
        function  wmvfn(me,src,evnt)
    %         co = get(src,'CurrentObject');
    %         if ismember(co,me.imhandles)
    %            set(src,'Pointer','fullcross')
    %         else
    %             set(src,'pointer','arrow')            
    %         end
        end
         %%%   
         function delete(me)
             if ishandle(me.fig) && ~strcmp(get(me.fig,'BeingDeleted'),'on')
                delete(me.fig) 
            
             end
         end
         %%%
         function destruc(me)
             for i = 1:length(me.sisobj)
                 me.sisobj(i).sisobj(me.sisobj(i).sisters==me) = [];
             end
             delete(me.extrafigs(ishandle(me.extrafigs)))
             delete(me);
         end
         %%%
         function meshButtonDown(me,src,evnt)

             meshhs = cat(1,me.meshes.ploth);
             actvm = find(any(meshhs == src,2));
             axh = get(src,'parent');
             ax = me.plotax([me.plotax.h]==axh);
             me.activemesh = actvm;
             me.activeax = ax;
             me.current_object = me.meshes(actvm);
 
             cp = get(axh,'currentpoint');

  
               me.bdownpt = ax.ax2vol(cp(1,1:2));
             if get(me.meshadjust,'value')==1
                 me.moving_mesh = true;
             end

         end
         %%%
         function meshButtonUp(me,src,evnt)
            if ~me.moving_mesh;
                return
            end
             ax = me.activeax;
  
%              cp = me.currentpoint;
            cp = get(ax.h,'currentpoint');
            msh = me.meshes(me.activemesh);

            me.buppt = ax.ax2vol(cp(1,1:2));
            dpos = me.buppt-me.bdownpt;
            msh.trirep = TriRep(msh.trirep.Triangulation,msh.trirep.X + repmat(dpos,size(msh.trirep.X(:,1))));
            me.meshes(me.activemesh) = msh;
            me.moving_mesh = false;
            me.plotmeshes;
         end
         %%%
         function figButtonUp(me,src,evnt)
             %%% Callback for button up
             if me.moving_mesh
                me.meshButtonUp(src,evnt);
             end
         end
        
%         %% Apply rotation matrix
%         function T = rtrmat(me,ax)
%              
%               T = ax.rtrmat;
%          
%  
%         end
        %%%%  Callback for the list menus    
        function listcallback(me,src,evnt)
            me.currlist = src;

            switch src

                case me.ptlist
                    fld = 'points';
                case me.meshlist
                    fld = 'meshes';
                case me.vollist
                    fld ='volumes';
                case me.transflist
                    fld ='transforms';
                case me.modulelist 
                    fld = 'modules';
   
                otherwise
                    error('unrecognized callback source')
            end
            indx = get(src,'value')-1;
            indx(indx==0) = [];
    %         me.current_object(1:end) = [];
       doubleclick = isequal(get(gcf,'selectiontype'),'open');
          if ~isempty(indx)
               switch fld
                   case 'modules'
                       if doubleclick
                           mod = me.modules{indx};
                              out= SubModules.(mod)(me); 
                              if isa(out,'SubModules.module')
                                   me.modobjects.(mod) = [me.modobjects.(mod),out];
                              end
                       end
                   otherwise

                           me.current_object= me.(fld)(indx);

                      if doubleclick %% Double click goes to object                            
                               me.gotoobj(); 
                      end
               end
            else
                me.current_object=me.info;
           end
            
            me.plotupdate();
            me.annotationUpdate(src,evnt);
        end
        %%%
        function ptboxupdate(me,src,evnt)
            i = me.pointbox==src;
            newval =str2double(get(src,'string'));
            cp = me.current_point;
            if ~isempty(newval)
                cp(i) = newval;
            end
            me.current_point = cp;

        end
            %%%
        function coordboxupdate(me,src,evnt)
            i = me.coordbox==src;
            newval =str2double(get(src,'string'));
            crT = me.current.transforms;
            cp = crT.tr(me.current_point);
            if ~isempty(newval)
                cp(i) = newval;
            end
            me.current_point = crT.itr(cp);

        end
     %%%            
        function modhelp(me,src,evnt)
           indx = get(me.modulelist,'value')-1;
           if indx == 0
               return
           else
               doc(['SubModules.',me.modules{indx}])
           end
           
        end

        %%%            
        function delobj(me,src,evnt,menuh)
           cro = me.current_object;
           fld = class(cro);
           switch fld
               case 'meshes'
                   lst = me.meshlist;
               case 'points'
                   lst = me.ptlist;
               case 'volumes'
                   lst = me.vollist;
               case 'transforms'
                   lst = me.transflist;
           end

           obji = find(ismember([me.(fld).objectid],[cro.objectid]));
           if all([cro.objectid] >0) % info field stays
               for k = length(obji):-1:1
                   if isa(cro,'plotobj')
                       plh = me.(fld)(obji(k)).ploth;
                       delete(plh(ishandle(plh)))
                   end
                   str = get(lst,'string');
                   str(obji(k)+1) = [];
                   set(lst,'string',str,'value',obji(k));
                   me.(fld)(obji(k)) = [];
               end
           end
        end
        %%%
        function updatelists(me)

                vis = {'-','o'};
                strfun = @(obj) [{upper(class(obj))},arrayfun(@(x)sprintf('%s   %s',vis{x.show+1},x.label),obj,'uniformoutput',false)];

                set(me.meshlist,'string',strfun(me.meshes));
                 set(me.ptlist,'string',strfun(me.points));
                 set(me.vollist,'string',strfun(me.volumes));
                 set(me.transflist,'string',{'TRANSFORMS',me.transforms.label});
                set(me.modulelist,'string',cat(1,{'MODULES'},me.modules));


        end
        %%% Callback for the annotation box
        function AnnoCallback(me,src,evnt)

            str = get(src,'string');
            cro = me.current_object;
            fld = class(cro);
            indx =ismember( [me.(fld).objectid] , [cro.objectid]);
            me.(fld)(indx).notes = str;        
        end

        %%%
        function gotoobj(me,src,evnt)

            cro = me.current_object;
            fld = class(cro);
            switch fld

                case 'meshes'
                    getpt  = @(x)mean(x.trirep.X,1);
                case 'points'
                    getpt  = @(x)x.coord;
                case 'volumes'
                    me.current.volumes = cro;
                    getpt = @(x)me.current_point;
                case 'info'               
                    return
                otherwise
                    error('unrecognized object')
            end
            if length(me.current_object)==1
                 indx = [me.(fld).objectid]==cro.objectid;
                 me.current_point = getpt(me.(fld)(indx));
            else
                pts = arrayfun(getpt,me.current_object,'uniformoutput',false);
                me.current_point = mean(cat(1,pts{:}));
            end
    %         me.plotupdate();

        end
        %%%
        function writemesh(me,src,evnt)

            cro = me.current_object;
            fld = class(cro);

            switch fld

                case 'meshes'
                    wrfun  = @(x,fn)arrayfun(@(x)writevtk(tri2vtk(x.trirep),[x.label,'.vtk']),x);
                    fn ='';
                case 'points'                
                    wrfun  = @(x,fn)dlmwrite([fn,'.txt'],cat(1,x.coord));
                    [fn,pth] =  uiputfile({'*.vtk','vtk';'*.txt','txt'},'What file to write?');
                    [~,fn,~] = fileparts(fn);
                otherwise
                    error('unrecognized object')
            end
%             indx = find(ismember([me.(fld).objectid],[cro.objectid]));

%             if length(indx)==1
                wrfun(cro,fullfile(pth,fn));
%             else
%                 for i = 1:length(indx)                
%                     wrfun(me.(fld)(indx),fullfile(pth,sprintf('%s%i',fn,i)));
%                 end
%             end
        end
        %%%%
        function hideobj(me,src,evnt)

            cro = me.current_object;
    %         fld = class(cro);
    %         oidx = find(ismember([me.(fld).objectid],[cro.objectid]));
                  for i = 1:length(cro)
                     cro(i).show = 1- cro(i).show;
                  end
            me.plotupdate();
        end  
            %%%%
        function renameobj(me,src,evnt)

            if ischar(src)
                str = [' for ',src];
            else
                str = '';
            end
            cro = me.current_object;
            fld = class(cro);
            indx = [me.(fld).objectid]==cro.objectid;

            newname =inputdlg(sprintf('Enter Label %s',str),'Enter Label',1,{me.(fld)(indx).label});
            me.(fld)(indx).label = newname{1};
            me.plotupdate();
        end 

        %%%%%

        function efigdel(me,src,evnt)
            ax = get(src,'children');
            axi =ismember([me.plotax.h],ax);
            me.plotax(axi) = [];
            me.extrafigs(me.extrafigs==src) = [];
        end
        %%%%%
        function updatePlotProps(me,src,evnt)
           co  = me.current_object;
           if ~isa(co,'plotobj'), return, end
           val = get(src,'string');
                    
           switch find(src==me.plprop);
               case 1
                    cln = str2num(val); 
                    if ~isempty(cln)
                        val = cln;
                    end
                    fld = 'plotcolor';
               case 2
                    fld = 'linestyle';
               case 3
                   fld = 'plotargs';
                   val = textscan(['x,',val],'%s','delimiter',',');
                   val = [val{:}];
                   val = val(2:end);
                   isn = cellfun(@(x)~isempty(str2num(x)),val); %#ok<*ST2NM>
                   val(isn) = cellfun(@(x)str2num(x),val(isn),'uniformoutput',false);
                   
           end
           for i = 1:length(co)
               co(i).(fld)= val;
           end            
        
        end
        
        function legendcheck(me)
           me.plotupdate(); 
        end
    end
end
 

