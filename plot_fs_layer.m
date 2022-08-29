function [c,cmap,unqi,names,vl] = plot_fs_layer(h,annot,show)


vl = [];
for k = 1:length(annot)
    
    [~,vlx,clt] = read_annotation(annot{k});
    vl = [vl;vlx];
end
if nargin < 3,
    show = true;
    mask = true(size(vl));
elseif isscalar(show)
    mask = 1;
else
    mask = show;
end

[ism,ismi]= ismember(vl.*mask,clt.table(:,end));
ismi(ismi==0)=1; %1st is usually "unknown".
[unq,~,unqi] = unique(ismi);
names = clt.struct_names(unq);

cmap = clt.table(unq,1:3)/255;

if isscalar(show) && ~show    
    c = [];
    return
end

% cmap =  [.5 .5 .5;cmap];
%set(h,'FaceVertexCdata',unqi ,'facevertexalpha',double(mask));
set(h,'FaceVertexCdata',cmap(unqi,:) ,'facevertexalpha',double(mask));
alpha('none')
caxis([.5 max(unqi)+.5])
colormap(cmap);

pos = get(gca,'position');
c = colorbar('west');
% cpos = get(c,'position');
% cpos(1) = .01;
% set(c,'position',cpos)

set(gca,'position',pos)


caxis([.5 max(unqi)+.5])
set(c,'ytick',1:size(clt.table(unq,:),1),'yticklabel',names)
  shading interp
%  shading flat

alpha flat