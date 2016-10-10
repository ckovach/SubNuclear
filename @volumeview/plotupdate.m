 function plotupdate(me,axs,updatesis)

 %%% Update Plots
 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%  me.cleanup();
 if nargin < 3
     updatesis = true;
 end

 if nargin < 2 || isempty(axs)
    axs = me.plotax; 
 elseif isnumeric(axs)
     axs = me.plotax(axs);
 end
 crp = me.current_point;
 getax = isopen(me.plotax);
 
 axs = axs(getax);
 
 if isempty(axs)
    warning('No axes to update...')
    return
 end
 
me.isplotting = true;
plh = cat(1,axs.ploth);
 for k = 1:length(me.volumes)
     if ~me.volumes(k).show
         set(plh(:,k),'visible','off')
         continue
     else
        for i = 1:numel(axs)
           axs(i).plotupdate(updatesis);

        end
     end


 end
 me.plotmeshes;
 arrayfun(@(a,b) set(a,'String',round(b)),me.pointbox,round(crp));
 arrayfun(@(a,b) set(a,'String',b),me.coordbox,round(me.current.transforms(1).tr(crp)));
 me.annotationUpdate();
   plhs = cat(1,me.meshes.ploth,me.points.ploth);
%    shown = cat(1,me.meshes.show,me.points.show)==1;
    able = ishandle(plhs);
    vis  =false(size(able));
    vis(able) =  strcmp(get(plhs(able),'visible'),'on');
   shown = any(vis,2);
   plhh = max(plhs'.*vis');
 if ishandle(me.legendh), delete(me.legendh), end
   
 if any(shown) && get(me.showlegendh,'value')==1;
     
         lbls = {me.meshes.label,me.points.label};
        l =  legend(me.plotax(1).h,plhh(shown),lbls(shown));
        pos = get(l,'position');
        pos(1:2) = [.5 .46];
        me.legendh = l;
        set(l,'units','normalized','location','none','position',pos ) 
  
%  else
%      legend off
 end
me.isplotting = false;
