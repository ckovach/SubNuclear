function trans = tpstransform(R1,R2,varargin)

% trans = tpstransform(R1,R2,arg3,arg4,...);
%
% Creates a transforms object for tps warping with control points R1 and
% R2. All additional arguments are passed to tpswarp.
%
% See also TPSWARP

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if nargin > 2 && isnumeric(varargin{1})
    reg = varargin{1};
    varargin(1) = [];
else
    reg = 0;
end

%%% Use SVD to handle degenerate cases
[u1,l1] = svd(cov(R1));
P1 = u1(:,diag(l1)>eps);
[u2,l2] = svd(cov(R2));
P2 = u2(:,diag(l2)>eps);

R1 = R1*P1;
R2 = R2*P2;

[~,trmat,nldf] = tpswarp(R1,R2,reg);

% trmat(size(R1,2)+1,size(R2,2)+1) = 1;
itrmat = trmat\eye(size(R1,2)+1);

pad = @(x)cat(2,x,ones(size(x,1),1));
iwarp= tpswarp(R2,R1); % Inverse warping. NB this is not a true inverse but only an approximation
                       % as it is not the case that iwarp(warp(x)) = x for
                       % all x. 
nlwarp = @(x)nldf(x*P1)*P2';
P1(end+1,end+1) = 1;                       
trans = transforms('trmat',P1*trmat*P2','nldfun',nlwarp,varargin{:},'type','TPS');


trans.inldfun = @(x)-pad(x*P2)*itrmat*P1(1:end-1,:)' + iwarp(x*P2);