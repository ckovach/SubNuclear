function [pts,transf,itransfull] = fitgrid(P,MN,PX,type,label,gridscale)

% pts = fit_grid(P,[m n],PX,type);
%
% Interpolates the points of an m x n grid based on the coordinates in matrix P.
%
% Inputs:
%
% P: Kx3 matrix
%
% [m n]: dimensions of the grid
%
% PX: Grid points corresponding to the rows of P 
%   If PX is not given or empty it defaults to
%       PX = [1 1; m 1; 1 n; m n]
%
% type: Uses linear interpolation if type = true and otherwise thin-plate spline
%       warping to project grid points into the space of P.
%      If type is not given or empty, it defaults to true when k < 5 and
%      otherwise false.
%
% Outputs:
%
% pts: m*n x 3 array of interpolated points.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2013

if nargin < 6
    gridscale = [1 1];
end

MN(3) = 1;
if nargin < 3 || isempty(PX)
    PX = MN([3 3; 1 3; 3 2; 1 2])*diag(gridscale);
end

if  nargin < 4 || isempty(type);
  type =  (size(PX,1) < 5)+1; 
end

if nargin  <5 || isempty(label)
   label = 'Grid %i'; 
end


if isa(P,'points')
    P = cat(1,P.coord);
end
[x, y] = meshgrid(1:MN(2),1:MN(1));

switch double(type)
    case 1
        P(:,4) = 1;
        PX(:,3) = 1;
        T = PX\P;
        transf = transforms('trmat',T,'label','To Grid');
        PX = PX(:,1:2);
        P = P(:,1:3);
    %     trfun = @(X) [X ones(size(X,1),1)]*T(:,1:end-1); 
    
   
    case 2
%        transf = transforms('label','To Grid');
%        tpsfun = tpswarp(PX,P,1e-9); 
       transf = tpstransform(PX*diag(gridscale),P,1e-9,'label','Grid2Vol');
      
    %    trfun = @(X) tpsfun(X);
%        transf.tr = tpsfun;
%        if nargout > 1
%            tpsifun = tpswarp(P,PX,1e-9); 
%           transf.itr = tpsifun;
%        end   
    case 3
        mp = mean(P);
        mpx = mean(PX);
        ZP = P-repmat(mp,size(P,1),1);
        ZPX = PX-repmat(mpx,size(PX,1),1);
        ZP(:,4) = 1;
        ZPX(:,3) = 1;
        T = ZPX\ZP;
        [u,l,v] = svd(T(1:2,1:3));
        l(1:2,1:2) = eye(2);
        T(1:2,1:3) = u*l*v';
        T = [eye(2,3);-mpx,1]*T*[eye(3,4);mp,1];
        transf = transforms('trmat',T,'label','To Grid');

        
end

pts = points(transf.tr([y(:) x(:) ]*diag(gridscale)));
for k = 1:length(pts)
    pts(k).label = sprintf(label,k);
end

%%% Compute an approximate inverse transformation
if nargout > 2
       PX0 = PX*diag(gridscale);
       PX0(:,3) = 0;
%        PX0 = [[y(:), x(:)]*diag(gridscale), zeros(size(x(:)))];
       PX1 = PX0;
       PX1(:,3) = mean(gridscale)*.1;
       nv = cross(transf.trmat(1,1:3),transf.trmat(2,1:3));
       nv = nv./norm(nv)*mean(gridscale)*.1;
        P0 = P;
       %        P0 = cat(1,pts.coord);
       P1 = P0 + repmat(nv,size(P0,1),1);
       itransfull = tpstransform(cat(1,P0,P1),cat(1,PX0,PX1),'label','Vol2Grid');
end
