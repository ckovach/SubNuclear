
function [p,facet] = slicetri(Tr,pt,nvec)

% function [p,facet] = slicetri(Tr,pt,nvec)
%
% Slice the TriRep object in the plane normal to nvec containing point pt.
%  
%  inputs:
%
%  Tr - is a mesh object describing a 2D manifold in a 3D euclidean space
%  pt - s point in the plane.
%  nvec - vector normal to the plane.
%
%  outputs
%
%  p - points at which the plane intersects facet edges, ordered according
%      to adjacency.
%  facet - indices into Tr.Triangulationn for the corresponding facets.
%  
%
%
% C Kovach 2013

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

nvec = nvec/norm(nvec); %make sure it's unit length



%%% If the plane intersects a facet, then one vertex must be on one side and
%%% the other two vertices on the other, which is reflected in the sign of the inner
%%% product with nvec.

X = Tr.X - repmat(pt,size(Tr.X,1),1);
prvs = sign(X*nvec');
prvs(prvs==0) = 1;
prvs = prvs(Tr.Triangulation);
fprvs = sum(prvs,2);  %%% if plane intersects, then this is sum is + or - 1, otherwise + or - 3
vppos = fprvs(:,[1 1 1])'.*prvs'; 
fint = abs(fprvs)~=3; %%% Facets that intersect

tri = Tr.Triangulation';

v0 = X(tri(vppos==-1),:);  %%% The single vertex on one side has value -1

v12 = X(tri(vppos==1),:);  %%% The two vertices on the other side each have value +1

v00 = zeros(size(v12));
v00(1:2:end,:) = v0;
v00(2:2:end,:) = v0;
d = v12-v00;

%%% These are the points at which the plane intersects the edges
xint = repmat(-v00*nvec'./(d*nvec'+eps),1,3).*d+v00;
% xint1 = repmat(-v0*nvec'./(d1*nvec'+eps),1,3).*d1+v0;
% xint2 = repmat(-v0*nvec'./(d2*nvec'+eps),1,3).*d2+v0;


% facets = Tr.Triangulation(fint,:);

%%% Now order points by neighbor.
fint(end+1) = 0;

ffint = find(fint(1:end-1));
ngh = Tr.neighbors(ffint);
edges = isnan(ngh);
ngh(isnan(ngh)) = length(fint);
nghin = fint(ngh)';
smng = sum(nghin);
if isequal(smng,0)
    p = [nan nan nan];
    facet = [];
    return
end

if any(smng == 3)
    try
    %     crpr = @(x,tr)cross(x(tr(:,2),:)-x(tr(:,1),:),x(tr(:,3),:)-x(tr(:,1),:));
    %compute area of a facet.
%     facearea = sqrt(sum(crpr(Tr.X,Tr.Triangulation).^2,2));
    dd = squareform(pdist(xint));
    [~,cri] = max(mean(dd));
    remap = zeros(size(xint,1),1);
    for k = 1:size(xint,1)
        remap(k) = cri;
        dd(cri,cri) = Inf;
        [~,mni] = min(dd(cri,:));
        dd(:,cri) = Inf;
        cri = mni;
    end
%     
%     [u,v] = svds(dd,3);
%     u = u(:,2:3)./repmat(sqrt(sum(u(:,2:3).^2,2)),1,2);
%     ang = atan2(u(:,1),u(:,2));
%     [srt,srti] = sort(ang);
% %     dd =squareform(pdist(xint));
% %     p = xint + repmat(pt,size(xint,1),1);
    p= xint(remap,:)+ repmat(pt,size(xint,1),1);
    facet = find(fint(1:end-1));

    catch
        p = xint;
         warning('Failed to reorder points. Points are unordered')
    end
else

% if any(smng==3)    
%     ff = find(smng==3);
%     fi = ff(ceil(rand*length(ff)));   
%     nhorig = ngh(ff,:);
% elseif any(smng~=-2)
%     ff = [];
%     fi = [];   
%     nhorig = [];
% 
% end
% 
% ci = 0;
% while any(smng==3) & ci < 500
% %         ng2 = Tr.neighbors(ngh(fi,nghin(:,fi))');
% %         smng2 = sum(fint(ng2),2) ;
%         nngi2 = ngh(fi,fint(ngh(fi,:)));
%         nngi2= nngi2(randperm(length(nngi2)));
%         fint(nngi2(1)) = 0;
%         
%         ffint = find(fint(1:end-1));
%         nghold= ngh;
%         ngh = Tr.neighbors(ffint);
%         ngh(isnan(ngh)) = length(fint);
%         edges = isnan(ngh);
%         nghin = fint(ngh)';
%         smng = sum(nghin);
%         if any(smng + sum(edges')==1) %#ok<*UDIM>
%              fint(nhorig') = fint(nhorig') | rand(size(nhorig'))<.0001;
% %             fint(ffint(ff),:) = rand(size(ff)*1;
%             fint(nngi2(1)) = 1;
%             ffint = find(fint(1:end-1));
%             ngh = Tr.neighbors(ffint);
%             ngh(isnan(ngh)) = length(fint);
%             edges = isnan(ngh);
%             nghin = fint(ngh)';
%             smng = sum(nghin);
%         end
%         if any(smng==3)
%             ff = find(smng==3);
%             fi = ff(ceil(rand*length(ff)));
%         end
% %         sum(smng==3)
%         ci = ci+1;
% end    
% % if mod(size(ngh,1),2)==1
% %     ngh(end+1,:) = [nan nan nan];
% % end

ngi = ngh';
ngx = ngi(nghin | edges');
ngi = reshape(ngx(1:size(ngh,1)*2),2,size(ngh,1),1);

if any(sum(edges,2)>0)
    indx = find(sum(edges,2)>0,1,'first');
    lsp = ngh(indx,edges(indx,:)); 
else
   indx = 1;
   lsp = ngi(1,1);
end

crp = ffint(indx);
p = nan(size(ngi,2),3);
if nargout > 1
    facet = nan(size(ngi,2),1);
end
for i = 1:size(ngi,2)
   
    indx = find(ffint==crp,1,'first');
    if isempty(indx)
        continue
    end
        
    nxpi = [1 2]*(ngi(:,indx)~=lsp);
    if nxpi == 3 || nxpi== 0
        nxpi = ceil(rand*2);
    end
        
    nxp = ngi(nxpi,ffint==crp);
    p(i,:) = xint(2*(indx-1)+nxpi,:)+pt;
    lsp = crp;
    crp = nxp;
    if nargout > 1
        facet(i) = crp;
    end
end
   
end


