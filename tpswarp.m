
function [warpfun,d,nlfun,D] = tpswarp(X1,X2,reg,dim)

% Warpfun = tpswarp(Xfrom,Xto, reg)
% Thin-plate spline warping. 
%      
%      Xfrom : row matrix of coordinates in the starting space.
%      Xto:   coordinates in the target space.
%      reg:  regularization term (default reg=0). For large values warping
%            approaches an affine transformation and for reg=0 it fits
%            the control points exactly.
%      dim: dimensionality of the warping space (default is the lower
%           dimensionality of X1 and X2)
%
%      Warpfun is a function handle where
%      Xto == Warpfun(Xfrom)   ( for reg = 0 );
% 
% [Warpfun,A] = tpswarp(Xfrom,Xto, reg)
%      Also returns the affine component of the transformation as matrix A.
%
% [Warpfun,A,Dfun] = tpswarp(Xfrom,Xto, reg)
%      Dfun is a function handle returning the non-linear displacement (warping minus the affine component) as an additional argument.
%


% C Kovach 2013

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if nargin < 3 || isempty(reg)
    reg = 0;
end

if isa(X1,'points')
    X1 = cat(1,X1.coord);
end
if isa(X2,'points')
    X2 = cat(1,X2.coord);
end

% squash = @(x)sqrt(sum(cat(3,x{:}),3));
% dxc = @(x,c,k) ;
% mdxc = @(x,c) squash(arrayfun(@(k)dxc(x,c,k),1:size(x,2)-1,'uniformoutput',false));
% tpsmat = @(x,c) mdxc(x,c).^2.*log(mdxc(x,c)+eps);
pad1 = @(x) cat(2,x,ones(size(x(:,1))));

if nargin < 4 || isempty(dim)
    dim = min(size(X1,2),size(X2,2));
end

%%% If the first input is of deficient rank, then project into the lower
%%% rank space
if rank(cov(X1))<min(size(X1))
    [u,l] = svd(cov(X1));
    u = u(:,diag(l)./l(1)>eps);
    u(end+1,end+1)=1;
    preProject = [eye(size(X1,2),size(X1,2)+1);-mean(X1) 1]*u;
else
    preProject = eye(size(X1,2)+1);
end
X1proj = pad1(X1)*preProject;

KX = @(x)tpsmat(x,X1proj(:,1:end-1),dim);

[Q,R] = qr(X1proj);
Q1 = Q(:,1:size(X1proj,2));
Q2 = Q(:,size(X1proj,2)+1:end);
 R = R(1:size(X1proj,2),:);
%R = R(1:rank(X1)+1,:);
TPScoef = Q2*(((Q2'*KX(X1proj(:,1:end-1))*Q2+reg*eye(size(Q2,2))))'\(Q2'*X2));
d = preProject*R^-1*Q1'*(X2-KX(X1proj(:,1:end-1))*TPScoef);
nlfun = @(x)KX(pad1(x)*preProject(:,1:end-1))*TPScoef;
warpfun = @(x) pad1(x)*d + nlfun(x);


 if nargout > 3
     D = @(x)mdxc(pad1(x)*preProject(:,1:end-1),X1*preProject(:,1:end-1));
 end

%%%%
function K = tpsmat(x,c,dim) 

% dim = size(x,2);
switch dim
    case {2,4}
        rbf = @(r) r.^(4-dim).*log(r+eps);
    otherwise
        rbf = @(r)r.^(4-dim);
end
% K = mdxc(x,c).^2.*log(mdxc(x,c)+eps); 
%         rbf = @(r) r.^2.*log(r+eps);
K = rbf(mdxc(x,c)); 
%%%%
function D = mdxc(x,c) 

d = 0;
% D = sqrt(squeeze(sum((repmat(x,[1, 1, size(c,1)])-repmat(permute(c,[3 2 1]),size(x,1),1)).^2,2)));
for k = 1:size(x,2)
    d = d+(repmat(x(:,k),1,length(c(:,k)))-repmat(c(:,k)',size(x,1),1)).^2;
%     D = D+d;
end
D = sqrt(d);

%  squash = @(x)sqrt(sum(cat(3,x{:}),3));
%  D = squash(arrayfun(@(k)dxc(x,c,k),1:size(x,2)-1,'uniformoutput',false));
  
