 
 function plotmeshes(me,varargin)
    
     % Plot meshes and points
     cols = prism;
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

     for k = 1:length(me.plotax)



          nv = me.plotax(k).normvec;
         for i = 1:length(me.meshes)
             vis = 'on';
             if ~me.meshes(i).show                     
                 slax = [0 0];
%                  set(me.meshes(i).ploth,'visible','off');
                  vis = 'off';
             else
                 sl = slicetri(me.meshes(i).trirep,me.current_point,nv);

                 slax = me.plotax(k).vol2ax(sl);
                 if isempty(slax)
                     vis = 'off';
                 end
             end
             if isempty(me.meshes(i).plotcolor)
                 me.meshes(i).plotcolor = cols(mod(i-1,length(cols))+1,:);
             end
             if isempty(me.meshes(i).linestyle)
                 me.meshes(i).linestyle = '-';
             end


             dslax = zscore(sum(diff(slax).^2,2));

             slax(dslax>10,:) = nan; %Large jumps are probably related to misordering

             set(me.meshes(i).ploth(k),'xdata',slax(:,1),'ydata',slax(:,2),'zdata',ones(size(slax(:,1)))*.5,'color',me.meshes(i).plotcolor,'linestyle',me.meshes(i).linestyle,'visible',vis,me.meshes(i).plotargs{:}) 
         end


         if ~isempty(me.points)                
             tol = 3; % how many pixels a point can be off the plain and still rendered
            for i = 1:length(me.points)
%             for i = find([me.points.show])
%                 try
               if me.points(i).show
                    pp = me.points(i).coord;
                    dp = pp-repmat(me.current_point,size(pp,1),1);
                    plp = abs(dp*nv') < tol;
%                     catch 
%                     warning('Deleted point %s with bad coordinate',me.points(i).label)
%                     delete(me.points(i))
%                     me.cleanup;
%                     break
%                 end
                    if any(plp) 
                         if isempty(me.points(i).plotcolor)
                             me.points(i).plotcolor = cols(mod(i-1,length(cols))+1,:);
                         end
                         if isempty(me.points(i).linestyle)
                             if size(me.points(i).coord,1) == 1

                                    me.points(i).linestyle = '+';
                             else
                                    me.points(i).linestyle = '-';

                             end
                         end

                        plx = me.plotax(k).vol2ax(pp);
                        plx = plx + 0./repmat(plp,1,2);
                        try
                            set(me.points(i).ploth(k),'xdata',plx(:,1),'ydata',plx(:,2),'visible','on','color',cols(mod(i-1,length(cols))+1,:),'zdata',.5,'markersize',14,...
                               'color',me.points(i).plotcolor,'linestyle',me.points(i).linestyle,me.points(i).plotargs{:})
                        catch err
                            warning(err.message)
                        end
                    else
                         try
                            set(me.points(i).ploth(k),'visible','off')
                         end
                    end
                else
                         try
                            set(me.points(i).ploth(k),'visible','off')
                         end
                end
            end
        end
     end


me.updatelists();

   