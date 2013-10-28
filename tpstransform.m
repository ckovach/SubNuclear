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

[~,trmat,nldf] = tpswarp(R1,R2,reg);

trmat(size(R1,2)+1,size(R2,2)+1) = 1;
itrmat = trmat\eye(size(R1,2)+1);

pad = @(x)cat(2,x,ones(size(x,1),1));
iwarp= tpswarp(R2,R1); % Inverse warping. NB this is not a true inverse but only an approximation
                       % as it is not the case that iwarp(warp(x)) = x for
                       % all x. 
                       
trans = transforms('trmat',trmat,'nldfun',nldf,varargin{:},'type','TPS');


trans.inldfun = @(x)-pad(x)*itrmat(:,1:end-1) + iwarp(x);