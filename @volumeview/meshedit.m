 function meshedit(me,mag, radius)
       
 %%% Edit meshes
 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if ~get(me.meshadjust,'value')
    warning('Press ''Adjust mesh'' radio buton to enable') 
    return
end

if nargin < 4
    radius =2;
end


crm = me.current.meshes(1);

Tr = me.current.meshes(1).trirep;

fn = Tr.faceNormals;
va = Tr.vertexAttachments;
vn = cellfun(@(x)mean(fn(x,:)),va,'uniformoutput',false);
vn = cellfun(@(x)x./norm(x),vn,'uniformoutput',false);
vn = cat(1,vn{:});

cp = me.current_point;

dv= @(x)sqrt(sum(x.^2,2));
D =Tr.X-repmat(cp,size(Tr.X,1),1);
d = dv(D)./(2*radius);
u = D./repmat(dv(D),1,3);        
fvec = u.*repmat(d.*exp(-(d).^2),1,3);

delta = repmat(sum(fvec.*vn,2),1,3).*u;
delta(isnan(delta)) = 0;
crm.trirep = TriRep(Tr.Triangulation,double(Tr.X + mag*delta));
me.plotupdate();
