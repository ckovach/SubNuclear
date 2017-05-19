
function [meshout,tr2tmpl] = coregister(Tr,templateTr,sphtr,affine)

% [meshresampled,tr2tmpl] = coregister(Tr,templateTr)
%
% Aligns and reorients two meshes, then projects vertices in templateTr
% onto Tr, creating a new mesh with the same connectivity as templateTr.
% 
%
%
% Tr - input mesh
% templateTr - template mesh
%
% meshresampld - retesselated mesh
% tr2tmpl affine transformation to template space;
%
%
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2013


if nargin < 3 | isempty(sphtr)
    addpath(fullfile(cd,'surfstat'))
    fprintf('Creating inflated mesh with call to freesurfer')
    surf.coord = templateTr.X';
    surf.tri= templateTr.Triangulation;
    SurfStatWriteSurf('template.surf',surf,'b')
    system('mris_inflate template.surf template_inflated.surf');
    sph = SurfStatReadSurf('rh.template_inflated.surf');
    
    sphtr = TriRep(double(sph.tri),sph.coord'); 
end
if nargin < 4
    affine = true; %% equate first and second moments if true
end

inc = sphtr.incenters;
[th,ph,r] = cart2sph(sphtr.X(:,1),sphtr.X(:,2),sphtr.X(:,3));
[sphc(:,1),sphc(:,2)] = cart2sph(inc(:,1),inc(:,2),inc(:,3));

%Cross product
crpr = @(x,tr)cross(x(tr(:,2),:)-x(tr(:,1),:),x(tr(:,3),:)-x(tr(:,1),:));
%Area of facet
facearea = @(tr)sqrt(sum(crpr(tr.X,tr.Triangulation).^2,2));
% Center of mass of surface
trcom = @(tr)facearea(tr)'*(tr.incenters./sum(facearea(tr)));
cnempt = @(C)~cellfun(@isempty,C);
plstruc = @(q)cellfun(@(x,k)trsurf(x,repmat(k,size(x.Triangulation,1),1)),q,num2cell(1:length(q)));
% %%% Thin plate spline
% dxc = @(x,c) (repmat(x,1,length(c))-repmat(c',length(x),1)).^2;
% mdxc = @(x,c) sqrt(dxc(x(:,1),c(:,1)));
% tpsmat = @(x,c) mdxc(x,c).^2.*log(mdxc(x,c)+eps);

%%% Find center of mass

mcent = trcom(Tr);
tmcent = trcom(templateTr);

% tmpTrcent = TriRep(templateTr.Triangulation,templateTr.X+repmat(mcent-tmcent,size(Tr.X,1),1));

%Now align second moments

cnt = @(x,w) x-repmat((w'*x)./sum(w),size(x,1),1);
wcov = @(x,w) cnt(x,w)'*(diag(w)./sum(w))*cnt(x,w);
trcov = @(tr)wcov(tr.incenters,facearea(tr));

TrCov = trcov(Tr);
tmpCov = trcov(templateTr);

shrink = @(x,a) a*(x-repmat(mean(x),size(x,1),1))+repmat(mean(x),size(x,1),1);

if affine
    % TRmat = (TrCov*sphCov^-1)^.5;
    %%% Rotate the principle components
    [umai,dmai] = svd(TrCov);
    [umni,dmni] = svd(tmpCov);
    umai = umai*diag(diag(sign(umai'*umni))); % Assuming they are roughly orinted the same
    TRmat = umni*(dmni^-1*dmai)^(1/2)*umai';
    template2tr = @(x) (x-repmat(tmcent,size(x,1),1))*TRmat + repmat(mcent,size(x,1),1) ;
    tr2templateMat = @(x) (x-repmat(mcent,size(x,1),1))*TRmat^-1 + repmat(mcent,size(x,1),1) ;
    template2trMat = [cat(1,TRmat,mcent-tmcent*TRmat),[0 0 0 1]'];
    tr2tmpl = template2trMat^-1;
    tr2tmpl = tr2tmpl(:,1:3);
    
%     trsphtr = TriRep(templateTr.Triangulation,shrink(template2tr(templateTr.X),.9));
else
    template2tr = @(x)x;
end
trsphtr = TriRep(templateTr.Triangulation,shrink(template2tr(templateTr.X),.9));
sphfn = trsphtr.faceNormals;
sphcn = trsphtr.incenters;

% for i =1:size(mnicn,1)
  
x1 = sphcn + 10*sphfn;  
x2 = sphcn - 10*sphfn;  

[in, on,pps] = intersects(x1,x2,Tr.Triangulation,Tr.X);  
in = in|on;
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
%%
%%%% Now interpolate vertex locations over spherical coordinates from the
%%%% inflated mesh.

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

X = vint;

meshout = TriRep(templateTr.Triangulation,X);


