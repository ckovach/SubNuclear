
function [meshout,sphTr] = isomesh(Tr,n,polyord,straightdim)

% [meshout,sphTr] = isomesh(Tr,n);
%
% Create a more homogenously spaced mesh by projecting an icosahedron onto
% the original mesh.
%
% Tr - input mesh
% n - mesh frequency
%
% meshout - retesselated mesh
% stphTr - transformed spherical mesh used to project vertices on to the
%            new mesh.
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2013

if nargin < 4
    straightdim = 2;  %%% Dimension along which to straighten
end

if nargin < 2
    n = 7; %frequency
end
if nargin < 3
    polyord = 8; %%% Order for polynomial straigtening
end

[icos(:,1),icos(:,2),icos(:,3),itr]= make_icosahedron(n,1,1);

micos = mean(icos);

sphtr = TriRep(itr,icos); 

% make sure that the face normals are directed out
inc = sphtr.incenters;
fn = sphtr.faceNormals;
fncov = sign(sum(fn.*inc,2));
itr(fncov<0,:) = itr(fncov<0,[1 3 2]);
sphtr = TriRep(itr,icos); 

%Cross product
crpr = @(x,tr)cross(x(tr(:,2),:)-x(tr(:,1),:),x(tr(:,3),:)-x(tr(:,1),:));
%Area of facet
facearea = @(tr)sqrt(sum(crpr(tr.X,tr.Triangulation).^2,2));
% Center of mass of surface
trcom = @(tr)facearea(tr)'*(tr.incenters./sum(facearea(tr)));
% cnempt = @(C)~cellfun(@isempty,C);
% plstruc = @(q)cellfun(@(x,k)trsurf(x,repmat(k,size(x.Triangulation,1),1)),q,num2cell(1:length(q)));
shrink = @(x,a) a*(x-repmat(mean(x),size(x,1),1))+repmat(mean(x),size(x,1),1);% %%% Thin plate spline
% dxc = @(x,c) (repmat(x,1,length(c))-repmat(c',length(x),1)).^2;
% mdxc = @(x,c) sqrt(dxc(x(:,1),c(:,1)));
% tpsmat = @(x,c) mdxc(x,c).^2.*log(mdxc(x,c)+eps);

%%% Find center of mass

mcent = trcom(Tr);

if straightdim ~=0
   %%% For irregularly shaped meshes, straighten along this dimension 
    Xc = Tr.incenters;
    xstr = Xc(:,straightdim);
    w = facearea(Tr);
    othdim =setdiff(1:3,straightdim);
    mkp = @(x)repmat(x,1,polyord+1).^repmat(0:polyord,length(x),1);
    Y = Xc(:,othdim);
    b = (diag(sqrt(w))*mkp(xstr))\(diag(sqrt(w))*Y);
    
    Xcent(:,straightdim) =Tr.X(:,straightdim)-mcent(straightdim);
    Xcent(:,othdim) = Tr.X(:,othdim)-mkp(Tr.X(:,straightdim))*b;
    Trcent = TriRep(Tr.Triangulation,Xcent);
else
    Trcent = TriRep(Tr.Triangulation,Tr.X-repmat(mcent,size(Tr.X,1),1));
end

%Now align second moments

cnt = @(x,w) x-repmat((w'*x)./sum(w),size(x,1),1);
wcov = @(x,w) cnt(x,w)'*(diag(w)./sum(w))*cnt(x,w);
trcov = @(tr)wcov(tr.incenters,facearea(tr));

TrCov = trcov(Trcent);
sphCov = trcov(sphtr);

TRmat = (TrCov*sphCov^-1)^.5;
% mai2t1amg = @(x) (x-repmat(maicent,size(x,1),1))*maiC^(-1/2)*t1C^(1/2) + repmat(t1cent,size(x,1),1) ;
% sph2mesh = @(x) x*TRmat + repmat(mcent,size(x,1),1) ;

in = false(size(Trcent.Triangulation,1),size(sphtr.Triangulation,1));
unmatched = sum(in)==0; % unmatched facets
for shrnk = 1:-.1:.1
     shrnk
    % progressively shrink the template until all facets are matched.
    trsphtr = TriRep(sphtr.Triangulation,shrink(sphtr.X*TRmat,shrnk));
    
    fns = trsphtr.faceNormals;
    sphfn(unmatched,:) = fns(unmatched,:);
    incs = trsphtr.incenters;
    sphcn(unmatched,:) =incs(unmatched,:);

    % for i =1:size(mnicn,1)

    x1 = sphcn(unmatched,:) + 10*sphfn(unmatched,:);  
    x2 = sphcn(unmatched,:) - 10*sphfn(unmatched,:);  

    [in(:,unmatched), on(:,unmatched),pps(:,unmatched,:)] = intersects(x1,x2,Trcent.Triangulation,Trcent.X);  
    in = in|on;
    unmatched = sum(in)==0; % Points that didn't find a dance partner
    if ~any(unmatched)
        break
    end
%         unmp = sphcn(unmatched,:);
%         dxc = @(x,c) (repmat(x,1,length(c))-repmat(c',length(x),1)).^2;
%         mdxc = @(x,c) sqrt(dxc(x(:,1),c(:,1))+dxc(x(:,2),c(:,2))+dxc(x(:,3),c(:,3)));
% 
%         dd = mdxc(Trcent.incenters,unmp);
%         [mn,mni] = min(dd);
%         in(mni+ (find(unmatched)-1)*size(in,1)) = true;
%     end
end
%%
for i = 1:size(sphcn,1)

  fac = in(:,i);
  pp = permute(pps(fac,i,:),[1 3 2]);

  dpp = sqrt(sum((pp-repmat(sphcn(i,:),size(pp,1),1)).^2,2));
  [mn,mni] = min(dpp);
  pp = pp(mni,:);      
  pproj(i,:) = pp;
  dproj(i) = mn;

end

%%%%

[th,ph,r] = cart2sph(icos(:,1),icos(:,2),icos(:,3));
[sphc(:,1),sphc(:,2)] = cart2sph(inc(:,1),inc(:,2),inc(:,3));

% xx = max(abs(ph))>max(abs(sphc(:,2)));

%%% Tile the azimuth and elevation since these are spherical coordinates
krmat = [0    0
         -2*pi 0 
          2*pi 0
         -pi pi
          pi pi
         -pi -pi
         pi -pi];
     
thphrep = [repmat(sphc,3,1);repmat(sphc*diag([1 -1]),4,1)] + kron(krmat,ones(size(sphc,1),1));
xyzrep = repmat(pproj,size(krmat,1),1);

dlt = DelaunayTri(thphrep);
for k = 1:3
  trint = TriScatteredInterp(dlt,xyzrep(:,k));
  vint(:,k) = trint(th,ph);
end

if straightdim>0
   X = vint;
   X(:,othdim) = X(:,othdim)+mkp(X(:,straightdim)+mean(xstr))*b;
   X(:,straightdim) =X(:,straightdim)+mcent(straightdim);
   
   if nargout >1
       sphX = trsphtr.X;
       sphX(:,othdim) = sphX(:,othdim)+mkp(sphX(:,straightdim)+mean(xstr))*b;
       sphX(:,straightdim) =sphX(:,straightdim)+mcent(straightdim);
       sphTr = TriRep(trsphtr.Triangulation,sphX);
   end
else
   X = vint + repmat(mcent,size(vint,1,1));
   if nargout >1
       sphTr = TriRep(trsphtr.Triangulation,trsphtr.X+repmat(mcent,size(trsphtr,X,1),1));
   end
end

meshout = TriRep(trsphtr.Triangulation,X);


