


 projection = 'ras_mni_aligned';

 hemi = 'lr';
 switch projection
     case 'ras_mni_aligned'
         %Compute a transform that maintains scanner coordinate scaling but rotates into
         %closer alignment with MNI space through a unitary transform.
         ras2mni = preop.transforms(1).trmat^-1*preop.transforms(3).trmat;
         [v,l,u] = svd(ras2mni(1:3,1:3));
         v(4,4)=1;u(4,4)=1;
         utrans = transforms('trmat',preop.transforms(1).trmat*v*u');
         voxtransform = @(x)utrans.tr(x);
     case 'ras'
         utrans = preop.transforms(1);
        voxtransform = @(x)preop.transforms(1).tr(x);
        
     case 'mni'
         utrans = preop.transforms(4);
        voxtransform = @(x)preop.transforms(3).tr(x);
     case 'vox'
         utrans = transforms('trmat',eye(4,3));
         voaxtransform = @(x)x;
     otherwise 
         error('Unrecognized projection %s',projection)
 end
  clear cxtrp
  
  addlink=false; 
%%
gethemi = ismember(hemis,hemi);
if ~exist('cxtrp','var') && exist('cxs','var')
    
    if length(cxs)>1 && sum(gethemi)>1
        if isa(cxs,'meshes')
            cxtrp = TriRep(cat(1,cxs(1).trirep.Triangulation,cxs(2).trirep.Triangulation+size(cxs(1).trirep.X,1)),voxtransform(cat(1,cxs(1).trirep.X,cxs(2).trirep.X))); %#ok<*DTRIREP>
            
        else
            cxtrp = TriRep(cat(1,cxs(1).tri,cxs(2).tri+size(cxs(1).vert,1)),voxtransform(cat(1,cxs.vert)));
        end
    else
        cxtrp = TriRep(cxs(gethemi).tri,voxtransform(cxs(gethemi).vert));
    end    
elseif ~exist('cxtrp','var')
    if ~isa(cx,'meshes')
        cxtrp = TriRep(cx.tri,vcx);
    else
        cxtrp = TriRep(cx.trirep.Triangulation,voxtransform(cx.trirep.X));
    end
end
one_liners

fig=figure;
h = trisurf(cxtrp);
shading flat, set(h,'FaceVertexCData',[1 .925 .925]*.75)
set(gcf,'color',[1 1 1]);

axis equal vis3d
light
material dull
 cameratoolbar


prp = preop_contact_locations;
lbl = {prp.label};
sf = @(x) ~cellfun(@isempty,strfind(lower(lbl),lower(x)));


%%% Pick which contacts to display
% geti = 1:length(preop_contact_locations);
% geti = (sf('grid') |sf('strip'))&~sf('sub')&~sf('parahip');
% geti = sf('right')&sf('cing')& abs(prx(:,1))'<10;
% geti =  sf('right');
%  geti = sf('strip')|sf('pole');
%       geti = (sf('sub')|sf('pole')|sf('parahip'))&~sf('depth')
%       geti = sf('depth')
%       geti = sf('right')
geti = true(size(lbl));
%        geti =  sf('grid') | sf('pole')| sf('strip')&~sf('sub') ;%& ~sf('subtemp') | sf('pole');
%      geti = sf('left')&~(sf('depth'))&~(sf('sub'))&~(sf('parahip'));
%  geti = sf('sub')|sf('pole')|sf('parahip')
%  geti = sf('posterior') & sf('temporal')
%      geti =  sf('grid') | sf('pole')| sf('strip')&~sf('sub') ;%& ~sf('subtemp') | sf('pole');
%     geti = sf('left')&~(sf('depth'))&~(sf('sub'))&~(sf('parahip'));
%    geti = sf('right')&(sf('grid') | sf('strip') | sf('interh') |sf('frontal')& ~sf('parahip') & ~sf('depth'));%& ~sf('sub')) ;
%    geti = sf('pole') | sf('sub') &  ~sf('depth') | sf('parahip');
%     geti =  sf('parahip') | sf('post') &  sf('sub') & sf('front');
%     geti=sf('grid') & (sf('_temporal') | sf('_frontal')); %| sf('strip');% sf('frontal') & ~sf('depth');
%      geti=sf('sub') | sf('pole') | sf('parahip'); %| sf('strip');% sf('frontal') & ~sf('depth');
%   geti=sf('strip') %| sf('strip');% sf('frontal') & ~sf('depth');
%    geti = sf('Heschl') | sf('insula')
% geti=sf('grid') | sf('strip') &~sf('sub')&~sf('parahipp');
%  geti = ~geti & ~sf('depth')
%geti=sf('superior temp') | sf('sub') ;
%geti = sf('insula');
% geti=sf('heschl') ;
%  geti=sf('frontal_grid') | sf('temporal_grid') | sf('pole') | sf('posterior') & sf('subfrontal') ;
%  geti=sf('occipital_grid') | sf('temporal_grid') | sf('pole') ;
% geti=sf('subocc') | sf('subtemp') | sf('pole') ;
%geti = sf('grid') & ~sf('subtemporal');
%geti = true(size(lbl));
% geti = sf('strip')&~sf('sub')&~sf('inter')&~sf('parah')&~sf('pole')|sf('grid');
%  geti= sf('grid') |sf('anterior_subfrontal')|sf('pole');
% geti = ~geti&~sf('depth')&~sf('interh');
% geti = sf('interh');
% geti=~sf('depth')&~sf('parahip');
% geti=sf('heschl');
%geti =  sf('depth')&diff([lozch.number,0])<0;
%geti2 = sf('grid')&sf('pole') &pabs(rx(:,1)')<=40;
% geti = sf('grid')&~sf('pole');% | geti2
%geti = sf('subtemporal') | sf('pole') | sf('subfrontal') | sf('parahippocampal') ;% | geti2
%geti = sf('depth')  ;% | geti2
% geti = sf('interhemispheric') | sf('subtemporal') | sf('parahip')  ;% | geti2

axis manual
prx = voxtransform(cat(1,prp.coord));
% hold on, plh =  mplot3(prx(geti,:),'k.','markersize',20,'linewidth',2);
hold on, plh =  scatter3(prx(geti,1),prx(geti,2),prx(geti,3),'k','fill');

%% After adusting rotation and lighting plot contacts
axis manual
shg
c = get(gca,'children');
delete(c(strcmp(get(c,'type'),'light')))
camlight headlight
lighting phong
axis off

projectview
% plh2 =  mplot3(pflat(geti,:),'ko','markersize',8,'linewidth',2);
% plh2 =  mplot3(pflat,'k.','markersize',20,'linewidth',2);
% plh2 =  mplot3(pflat,'k.');
plh2 =  scatter3(pflat(:,1),pflat(:,2),pflat(:,3),'k','fill');
%hold on, plh2 = mplot3(pflat(geti,:),'ko','markersize',15,'linewidth',3)
set(gca,'clipping','off')
 %% Save results
 saveas(gcf,fullfile(ddir,'brain.fig'))
 
axis on
grid on
print(gcf,'-dpng','-r400',fullfile(ddir,'brain_layer_w_contacts.png')) 
axis off

set([plh,plh2],'visible','off')

print(gcf,'-dpng','-r400',fullfile(ddir,'brain_layer.png'))

set(h,'visible','off')
set(plh2,'visible','on')

if addlink
%     txh = arrayfun(@(x,y,z,c)text(x,y,z,strtok(c.label,':'),'horizontalalignment','center'),pflat(:,1),pflat(:,2),pflat(:,3),prp(geti)','uniformoutput',false);
     txh = arrayfun(@(x,y,z,c)text(x,y,z,'  ','horizontalalignment','center'),pflat(:,1),pflat(:,2),pflat(:,3),prp(geti)','uniformoutput',false);
    txh=[txh{:}];
    baseurl = 'https://saccade.neurosurgery.uiowa.edu/labwiki/index.php/';
    hl = strcat(baseurl,arrayfun(@(x)sprintf('Subject/%s/Contact_%s',sid,strtok(x.label,':')),prp(geti),'uniformoutput',false));
    pdflink(fig,txh,hl,fullfile(ddir,'contact_layer.pdf'),'-r400');
else
    print(gcf,'-depsc','-r400','-painters',fullfile(ddir,'contact_layer.eps'))
end

set(h,'visible','on')
show = true;
%% Save freesurfer surface parcellations
shg
set([plh,plh2],'visible','off')

%  annot = fullfile(freesurfer_directory,'label','lh.BA.annot')
annot = {};
ghemis = hemis(gethemi);
for k = 1:length(ghemis)
    hmi = ghemis(k);
   annot{k} = fullfile(freesurfer_directory,'label',[hmi,'h.aparc.annot'])
end 
% annot = fullfile(freesurfer_directory,'label',[hemi,'h.aparc.DKTatlas40.annot']);
shading interp
c = plot_fs_layer(h,annot,show);
% set(h,'colorshading','interp')
set(c,'visible','off')
print(gcf,'-dpng','-r400',fullfile(ddir,'DKT_atlas_layer.png'));
set(c,'visible','on')
set(h,'visible','off')
print(gcf,'-depsc','-painters','-r400',fullfile(ddir,'DKT_atlas_layer_key.eps'));
delete(c)

set(h,'visible','on')
for k = 1:length(ghemis)
    hmi = ghemis(k);
     annot{k} = fullfile(freesurfer_directory,'label',[hmi,'h.aparc.a2009s.annot']);

end 
c = plot_fs_layer(h,annot,show);
set(c,'visible','off')
print(gcf,'-dpng','-r400',fullfile(ddir,'Destrieux_atlas_layer.png'));
set(c,'visible','on')
set(h,'visible','off')
print(gcf,'-depsc','-painters','-r400',fullfile(ddir,'Destrieux_atlas_layer_key.eps'));
delete(c)
set(h,'visible','on')

%% Project a selected point on the surface into the subject's image space
% When a point isset selected on the surface, the preop volumeview object is
% updated with the coordinates. Right clicking saves the selected point.
set(fig,'WindowButtonDownFcn',@(obj,hit)surf2subnuclear(obj,hit,preop,@(x)utrans.itr(x))); %Object callbacks will keep the surface and volumeview objects in sync
set(h,'ButtonDownFcn',@(obj,hit)surf2subnuclear(obj,hit,preop,@(x)utrans.itr(x)));
figure(preop.fig)
figure(fig)
hlpd = helpdlg('Pick a point on the surface');
set(hlpd,'units','normalized','position',[0.17 0.60 0.14 0.077]);
%% Projected images
shg
set(h,'visible','off')
colorbar off
[x2m,map] = rgb2ind(x2,1024);
caxis([0 1023])
decim = 4;
sh= surf(pproj2(1:decim:end,1:decim:end,1),pproj2(1:decim:end,1:decim:end,2),pproj2(1:decim:end,1:decim:end,3),double(x2m(1:decim:end,1:decim:end,:)));
set(sh,'cdatamapping','scaled')
shading flat
material dull
colormap(map)

print(gcf,'-dpng','-r400',fullfile(ddir,'implantation_image_layer.png'));

set(sh,'visible','off')
sh2= surf(pproj2(1:decim:end,1:decim:end,1),pproj2(1:decim:end,1:decim:end,2),pproj2(1:decim:end,1:decim:end,3),double(x2m(1:decim:end,1:decim:end,:))+0./inp3(1:decim:end,1:decim:end,:));
shading flat
material dull
print(gcf,'-dpng','-r400',fullfile(ddir,'implantation__cropped_image_layer.png'));%,'Alpha',inp3(1:decim:end,1:decim:end));
im = imread(fullfile(ddir,'implantation__cropped_image_layer.png'));
imwrite(im,fullfile(ddir,'implantation__cropped_image_layer.png'),'Alpha',1-double(all(im==255,3)));
% set(sh,'visible','off')
% sh2= surf(pproj2(:,:,1),pproj2(:,:,2),pproj2(:,:,3),double(x2)./255+0./inp3(:,:,[1 1 1]));
% shading flat
% material dull
% print(gcf,'-dpng','-r400',fullfile(ddir,'implantation__cropped_image_layer.png'));
% set(sh,'alphadatamapping','direct')


[x1m,map] = rgb2ind(x1(1:5:end,1:5:end,:),1024);

set([sh,sh2],'visible','off')
sh3= surf(pproj(:,:,1),pproj(:,:,2),pproj(:,:,3),double(x1m)+0./inp);
shading flat
material dull
print(gcf,'-dpng','-r400',fullfile(ddir,'implantation__cropped_image1_layer.png'));%,'Alpha',inp3(1:decim:end,1:decim:end));
im = imread(fullfile(ddir,'implantation__cropped_image1_layer.png'));
imwrite(im,fullfile(ddir,'implantation__cropped_image1_layer.png'),'Alpha',1-double(all(im==255,3)));

%% Plot superior temporal plane
hemi='l';
lats = {'right','left'};
lats = lats(1 + mod([1 2] + find(hemis=='r'),2));
%annot = fullfile(freesurfer_directory,'label',[hemi,'h.aparc.DKTatlas40.annot']);

annot = fullfile(freesurfer_directory,'label',[hemi,'h.aparc.annot']);

if length(cxs)>1
    hi = find(hemi==hemis);
    cxtrp = TriRep(cxs(hi).tri,voxtransform(cat(1,cxs(hi).vert)));
    vcxtrp = cxtrp;
end

[vi,vl,clt] = read_annotation(annot);

getstr = ismember(clt.struct_names,{'superiortemporal','transversetemporal','temporalpole'});
hgi = ismember(clt.struct_names,{'transversetemporal'});
show = ismember(vl,clt.table(getstr,5));
showhg = ismember(vl,clt.table(hgi,5));
hgface = any(showhg(vcxtrp.Triangulation),2);
% sttrp = TriRep(vcxtrp.Triangulation(any(show(vcxtrp.Triangulation),2),:),preop.transforms(3).tr(vcxtrp.X));
 sttrp = TriRep(cxtrp.Triangulation(any(show(vcxtrp.Triangulation),2),:),cxtrp.X);
% set(h,'facevertexalphadata',double(show))

fig=figure
h = trisurf(sttrp);
set(h,'facevertexalphadata',double(show))
shading flat, 
cdata = ones(size(sttrp.Triangulation,1),1)*[1 .925 .925]*.75;
%cdata(hgface,:) = cdata(hgface,:)*.7;

set(h,'FaceVertexCData',cdata)
set(gcf,'color',[1 1 1]);

axis equal vis3d
camlight
material dull
cameratoolbar
axis off
hold on

% lbl = {prp.label};
% gprp = prp(geti);
% sf = @(x) ~cellfun(@isempty,strfind(lbl,x));


%%% Pick which contacts to display
% sf = @(x) ~cellfun(@isempty,strfind(lbl,x));
geti = (sf('heschl')|sf('planum')| (sf('temporal')&sf('superior'))) & sf(lats{hemi==hemis}) ;
geti = (sf('heschl')|sf('planum')| (sf('temporal')&sf('superior')))  ;
% geti = sf(lats{hi})&(sf('insula') |sf('Heschl') |sf('temporal')&sf('superior'));% & sf('depth') & sf(lats{gethemi});
%  geti = sf(lats{hi})&(sf('insula') |sf('Heschl') |sf('temporal')&sf('superior')|(sf('cing')&sf('supra')))&~sf('oblique');% & sf('depth') & sf(lats{gethemi});
%  geti = cnum>1000&sf(lats{hi})&(sf('insula') |sf('Heschl') |sf('temporal')&sf('superior')|(sf('cing')&sf('supra')))&~sf('oblique');% & sf('depth') & sf(lats{gethemi});
%geti = (sf('insula') |sf('Heschl')|(sf('post')&sf('cing')) ) & sf('depth');
% geti = sf('grid')&~sf('pole');% | geti2
%geti = sf('subtemporal') | sf('pole') | sf('subfrontal') | sf('parahippocampal') ;% | geti2
% geti = sf('heschl') ;%| sf('hippocamp') & sf('depth')  ;% | geti2
% geti = sf('heschl') | sf('planum') ;%| sf('hippocamp') & sf('depth')  ;% | geti2

% hold on, plh =  mplot3(preop.transforms(3).tr(prx(geti,:)),'ko','markersize',10,'linewidth',2);
hold on, plh =  mplot3(prx(geti,:),'ko','markersize',8,'linewidth',2);

%%  Sagittal, coronal and axial views for all depth electrodes

mkdir(fullfile(ddir,'Depth_indiv'))
%geti = sf('depth')  ;% | geti2
% geti = 1:length(preop.points);
for k = 1:length(preop.points)
    preop.points(k).show = false;
    postop.points(k).show = false;    
end
h = helpdlg('Select points to plot in preop and postop images and press OK');
uiwait(h);

gpop = postop.current.points;
gprp = preop.current.points;

% glbl = lbl(geti);
glbl = {gprp.label};


subdir = fullfile(ddir,'Depth_indiv');

% gprp = preop.points(geti);
% gpop = postop.points(geti);


%for k = 1%:length(gprp)
fig = figure('renderer','painters','PaperSize',[11 8],'paperposition',[0 0  11 8])

set(fig,'color','k','inverthardcopy','off')
axhs = arrayfun(@(k)subplot(2,3,k),1:6,'uniformoutput',false); 
axhs =[axhs{:}];
ax2= axes; axtt = title(ax2,'','color','w','interpreter','none');
axis(ax2,'off')

axperm = [ 2 3 1 4
           1 3 2 4
           1 2 3 4];
 axtr = utrans.trmat*diag([-1 1  1 1]);
clear axpr axpo
for axi = 1:3
    arrayfun(@(h)cla(h),axhs(axi+([0,3])));
    axpr(axi) = preop.addaxis(axhs(axi),false);
    axpr(axi).showcrossHere=false;
    axpr(axi).showcrossThere=false;
    axis image
%     axpr(axi).trmat = preop.plotax(axi).trmat;
    axpr(axi).trmat = axtr(:,axperm(axi,:)); 
    
    pos = get(axpr(axi).h,'position');
    pos(1:2) = pos(1:2)-[0  .08];
    pos(3:4) = [.30 .4];
    set(axpr(axi).h,'position',pos);
    
    axpo(axi) = postop.addaxis(axhs(axi+3));
    axpo(axi).showcrossHere=false;
    axpo(axi).showcrossThere=false;
    axis image
%     axpo(axi).trmat= postop.plotax(axi).trmat;
    axpo(axi).trmat = axpr(axi).trmat; 
    pos = get(axpo(axi).h,'position');
    pos(1:2) = pos(1:2)-[0 .08];
    pos(3:4) = [.3 .4];
    set(axpo(axi).h,'position',pos);
    
    tpr(axi) = title(axpr(axi).h,'');
    tpo(axi) = title(axpo(axi).h,'');
    
%     
%     preop.current_point = gprp(k).coord;
%     postop.current_point = gpop(k).coord;
%     xlim(xlim)
%     ylim(ylim)
    %         drawnow
    %            axpo(axi).showcrossHere=true;
    %              axpr(axi).showcrossHere=true;
    %
    
  axpr(axi).scalebar.length=10;
  axpo(axi).scalebar.length=10;
  
end
preop.current_point = gprp(1).coord;
postop.current_point = gpop(1).coord;
%preop.resetaxis();
hd = [axpr.handles,axpo.handles];
hdh = [hd.control,hd.rotate,hd.transpose];
set(hdh,'visible','off')
axis([axpr.h,axpo.h],'image')

if ~exist('post_implant_label','var')
    post_implant_label='POST SURGICAL MRI';
end
    
%%
set(fig,'PaperSize',.75*[11 8],'paperposition',.75*[0 0  11 8])
drawnow
sufx = '';
shg
% sufx = '_preexp';
for k = 1:length(gprp)
%     for k = find((sf2('heschl')))
    
    gprp(k).show = true;
    gpop(k).show = true;
    preop.current_point = gprp(k).coord;
    postop.current_point = gpop(k).coord;
    
%     for axi =1:3
%         axpr(axi).showcrossHere=true;
%         axpo(axi).showcrossHere=true;
%     end
    gprp(k).linestyle = '+';
    gprp(k).plotargs = {'linewidth',.75,'markersize',8};
    gprp(k).plotcolor = [1 0 0];
%     gprp(k).marker= 'o';
        
    gpop(k).linestyle = '+';
    gpop(k).plotcolor = [1 0 0];
%     gpop(k).marker= 'o';  
    gpop(k).plotargs = {'linewidth',.75,'markersize',8};
    
         
    preop.plotupdate();
    postop.plotupdate();
   
%        delete([axpr,axpo])
    trs = preop.transforms([1 ]);
    prcoords = arrayfun(@(x)x.tr(gprp(k).coord),trs,'uniformoutput',false);
    arrayfun(@(x,y,z)set(x,'string',sprintf('%s: %0.0f mm',z,y)),tpr,prcoords{1},'XYZ')

    yl = ylabel(axpr(1).h,sprintf('PRE SURGICAL MRI'));%      \n-------------------------\nVoxel:     %0.0f, %0.0f, %0.0f\n\nRAS(mm):  %0.0f, %0.0f, %0.0f\n\nMNI(mm):  %0.0f, %0.0f, %0.0f',gprp(k).coord,trs(1).tr(gprp(k).coord),trs(3).tr(gprp(k).coord)),...
%                 'units','normalized','position',[-.1    0.5         0],'horizontalalignment','right');
    trs = postop.transforms([1 ]);
    pocoords = arrayfun(@(x)x.tr(gpop(k).coord),trs,'uniformoutput',false);
    arrayfun(@(x,y,z)set(x,'string',sprintf('%s: %0.0f mm',z,y)),tpo,pocoords{1},'XYZ')
   
   
    yl(2) = ylabel(axpo(1).h,post_implant_label);%\n-------------------------\nVoxel:     %0.0f, %0.0f, %0.0f\n\nRAS(mm):  %0.0f, %0.0f, %0.0f\n\nMNI(mm):  %0.0f, %0.0f, %0.0f',gpop(k).coord,trs(1).tr(gpop(k).coord),trs(3).tr(gpop(k).coord)),...
%                 'units','normalized','position',[ -.1   0.5         0],'horizontalalignment','right');
    axis([axpr(1).h,axpo(1).h],'on')
    set([axpr(1).h,axpo(1).h],'xtick',[],'ytick',[],'color','none')
    set(axtt,'string',glbl{k},'color','w')
    set([yl,tpr,tpo],'color','w')
    
    %%
  %  set(yl,'rotation',0)
    %th = title(axpr(2).h,regexprep(glbl{k},'_',' '),'interpreter','none','units','normalized','position',[.5 1 0])
    fn = fullfile(subdir,regexprep(glbl{k},'[:\s/]+','_'));
   % print(fig,'-dpng','-r400',[fn,'.png']);    
     print(fig,'-depsc','-r400',[fn,'.eps']);
%   print(fig,'-dpdf','-r200',[fn,'.pdf']);
%     [err,out]=system(sprintf('gs  -sDEVICE=png16m -dEPSCrop -r300 -o %s.png %s.eps',fn,fn));
    
   fns{k} = fn;
    gprp(k).show = false;
    gpop(k).show = false;

end
%%
fh = figure;
set(fh,'color','k','inverthardcopy','off')
axis off
[err,un{1}] = system('whoami');
if ~exist('hemi','var')
    hemi='';
end
th = text(.5,.7,sprintf('Patient %s%s contact locations',sid,upper(hemi)),'color','w','horizontalalignment','center','fontsize',24);
th2 = text(.5,.4,sprintf('Created:  %s\n\nby user:  %s',datestr(now),un{~err}),'color','w','horizontalalignment','center','fontsize',18);
fn = [fullfile(subdir,'title')];
% print(fh,'-dpdf','-r200',fn);
print(fh,'-depsc','-r400',fn);
 
% gt = gtext(sprintf('Patient %s%s contact locations',sid,upper(hemi)),'color','w','horizontalalignment','center','fontsize',30)
% gt2 = gtext(sprintf('\n\n\n\nCreated:  %s\n\nby user:  %s',datestr(now),un{~err}),'color','w','horizontalalignment','center','fontsize',16)

% infiles = sprintf(' %s.eps ',fn,fns{:});
infiles = sprintf(' %s.eps ',fn,fns{:});
infiles = regexprep(infiles,'\(','\\(');
infiles = regexprep(infiles,'\)','\\)');
com1 = sprintf('gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s%s%s_all_contacts%s.pdf %s',...
                           subdir,filesep,sid,sufx,infiles);
%  [err,out]=system(com1);
outfn = sprintf('%s%s%s_all_contacts%s.pdf',subdir,filesep,sid,sufx);
com2=sprintf('pdfcrop %s %s',outfn,outfn);
%  [err2,out2]=system(com2);

fid = fopen('temp.sh','w');
fprintf(fid,'LD_LIBRARY_PATH=\n%s\n\n%s',com1,com2);
fclose(fid);
system('chmod 755 temp.sh');
[err,out] = system('xterm -e "./temp.sh; read"');
%  [err,out]=system(sprintf('gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s%s%s_depth_contacts.pdf %s',...
%                            subdir,filesep,sid,infiles));
%                        
                       
%  [err,out]=system(sprintf('cp %s%s%s_depth_contacts.pdf ~/hbrl2/HBRL_Upload/Anatomical\ reconstruction/%s',...
%                            subdir,filesep,sid,sid));

%%
figure(fig)
hold(gca,'on')
if ~exist('plh','var')
    plh = [];
end

mtext = @(x,varargin)text(x(:,1),x(:,2),x(:,3),varargin{:});
lbl2 = arrayfun(@(x)sprintf('%i: %s %i',x.contact,x.label,x.number),blkdat.lozchannels,'uniformoutput',false)
sf2 = @(x) ~cellfun(@isempty,strfind(lower(lbl2),lower(x)));
figure(fig)
if any(sf2(postop.current.points(1).label))
    chan = blkdat.lozchannels(sf2(postop.current.points(1).label));
    chans = blkdat.lozchannels(sf2(chan.label));
else 
    chans=struct('contact', [1 length(postop.current.points)]);
    chan.label = postop.current.points(1).label;
end
plh(end+1) = mplot3(voxtransform(po2pr(postop.current_point)),'ro','markersize',10,'markerfacecolor','r');
plh(end+1) = mtext(voxtransform(po2pr(postop.current_point)),sprintf('%s (%i-%i)',chan.label,minmax([chans.contact])),'FontSize',8,'verticalalignment','bottom');
