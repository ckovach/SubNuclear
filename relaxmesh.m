function newmesh = relaxmesh(orig_mesh,template_mesh,curve_threshold)

%% relax the mesh so that vertex spacing approximates that of the template
niter = 10;
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

nrep = 5;
% dE = @(d,d0) d'*(1-d0./sqrt(sum(d.^2,2)));
% d2E = @(
% lnfy = @(x)x(:);
% trln = @(x)lnfy(x');
if nargin <3
    curve_threshold = .5; %throw out any vertices whos curvature exceeds some value
end
one_liners
newmesh = orig_mesh;
edgi = template_mesh.edges;
edgd0 = sqrt(sum((template_mesh.X(edgi(:,1),:)-template_mesh.X(edgi(:,2),:)).^2,2));
   fh = figure;
   ax1 = subplot(1,2,1);
   ax2 = subplot(1,2,2)
   E = []
fa = facearea(newmesh);
%  constraint = 'fn';
 constraint = 'vn';
%  method = 'grad';
%  method = 'jitter';
 method = 'grad'
 njitt = 100;
%  fr = getframe
%  fr=fr([]);
 E = [];
 
%  TRI = newmesh.Triangulation;
 TRI = newmesh.Triangulation;
 TRI(fa==0,:) = []; % Ignore faces with zero area
 newmesh = TriRep(TRI, newmesh.X);
  z2nan = @(x) x + 0./x; %convert
m2plot3 =   @(x,y,varargin)plot3(cat(2,x(:,1),y(:,1))',cat(2,x(:,2),y(:,2))',cat(2,x(:,3),y(:,3))',varargin{:})

jt = jet;
coli =  @(x)round((x-min(x))./diff(minmax(x))*(size(jt,1)-1)+1);

for kk = 1:nrep
 
    for i = 1:niter
    %     inc = orig_mesh.incenters;
        edgd = sqrt(sum((newmesh.X(edgi(:,1),:)-newmesh.X(edgi(:,2),:)).^2,2));
        E(end+1) = sum((edgd-edgd0).^2);
        axes(ax1)
        plot(E)
        xlim([0 nrep*niter]);
        i

        axes(ax2)
        pp = m2plot3(newmesh.X(edgi(:,1),:),newmesh.X(edgi(:,2),:));
        arrayfun(@(pp,col)set(pp,'color',jt(col,:)),pp,coli(edgd))
%          sm = arrayfun(@(xi,x) z2nan(ismember((1:size(newmesh.X,1))',xi{1})*(x+eps)),newmesh.edgeAttachments(edgi),(edgd-edgd0).^2,'uniformoutput',0);
%         edgE = nanmean([sm{:}],2);
%         trmesh(newmesh,edgE)
        colorbar
        caxis(minmax(edgd))
%         fr(end+1) = getframe(gcf);
        Xnew = newmesh.X;
        rv = randperm(size(newmesh.X,1));
    %     fns = newmesh.faceNormals;
    %     va = newmesh.vertexAttachments;
    %     vns= cellfun(@(x)mean(fns(x,:)),va,'uniformoutput',0);
    %     vns = cat(1,vns{:});

        for ki = 1:size(newmesh.X,1)

            k = rv(ki);

            v = newmesh.X(k,:);
            vorig = template_mesh.X(k,:);
            vtr = newmesh.vertexAttachments(k);
%             tri = TRI(vtr{1}',:);
            tri = newmesh.Triangulation(vtr{1}',:);
            vtr{1}(all(tri~=k,2)) = [];
            tri(all(tri~=k,2),:) = [];
            if length(vtr{1})<3
                continue
            end
            inc = newmesh.incenters(vtr{1}');
            fn = newmesh.faceNormals(vtr{1}');
            vn = mean(fn);
            vn = vn./norm(vn);
            
            [vatti,~,unqi] = unique(newmesh.Triangulation(vtr{1},:)');
%             [vatti,~,unqi] = unique(TRI(vtr{1},:)');
            unqi(unqi==find(vatti==k)) = [];
            unqi(unqi>find(vatti==k))=unqi(unqi>find(vatti==k))-1;
            unqi = reshape(unqi,2,size(tri,1))';
            vatti(vatti==k) = [];

            vatt = newmesh.X(vatti,:);
            vattorig = template_mesh.X(vatti,:);

            dvatt = @(v) repmat(v,size(vatt,1),1)-vatt; 
            dv1 = vatt(unqi(:,1),:)-repmat(v,size(unqi,1),1);
            dv2 = vatt(unqi(:,2),:)-repmat(v,size(unqi,1),1);
            dv3 = cross(dv1,dv2);
            dv3 = diag(sqrt(sum(dv3.^2,2)).^-1)*dv3;
            catdv = reshape(cat(1,dv1',dv2',dv3'),3,size(dv1,1)*3);
            % Project into the plane of the vertex normal
            catdvprj = catdv - repmat(vn',1,size(catdv,2))*diag(vn*catdv);
            catdvprj(:,3:3:end) = repmat(vn',1,size(dv1,1));
            dv = sparseblock(catdvprj,3);

            d0 = sqrt(sum((repmat(vorig,size(vattorig,1),1)-vattorig).^2,2)); 
            rr = @(v,d0) d0./sqrt(sum(dvatt(v).^2,2));       
            dE = @(v,d0) 2*dvatt(v)'*(1-rr(v,d0));
            %%% abs to make sure its convex
            d2E = @(v,d0) 2*sum(abs(1-rr(v,d0)))*eye(3) + 2*dvatt(v)'*diag(rr(v,d0).^3./d0.^2)*dvatt(v) ;

            % Constrain the vertices to move only in the facet planes.
    %         dEprjs = repmat(dE(v,d0),1,size(dv1,1)) - fn'*diag(fn*dE(v,d0));
    %         fprintf('\n')
            switch method
                case 'grad'
                    grE = -(d2E(v,d0)+0*eye(3))^-1*dE(v,d0);
                case 'jitter'
                    DE = @(x,y,z) sum((sqrt(sum(dvatt(v+[x,y,z]).^2,2))-d0).^2)-sum((sqrt(sum(dvatt(v).^2,2))-d0).^2);
                    rd =  (d2E(v,d0))^-1*trnd(1,3,njitt) ;
                    des = arrayfun(DE,rd(1,:),rd(2,:),rd(3,:));
                case 'grjitt'
                    grE = -(d2E(v,d0)+0*eye(3))^-1*dE(v,d0);
                    DE = @(x,y,z) sum((sqrt(sum(dvatt(v+[x,y,z]).^2,2))-d0).^2)-sum((sqrt(sum(dvatt(v).^2,2))-d0).^2);
    %                 DE1 = @(x,y,z) sum((sqrt(sum(dvatt(v+[x,y,z]).^2,2))-d0).^2)-sum((sqrt(sum(dvatt(v).^2,2))-d0).^2);
    %                 gr2 = @(x,y,z) sum(dvatt(v+[x,y,z])*vn'./sqrt(sum(dvatt(v+[x,y,z]).^2,2))); %try to minimize 2nd derivative too
    %                 DE = @(x,y,z)DE1(x,y,z) - .1*gr2(x,y,z);
                    rd =  (d2E(v,d0))^-1*trnd(1,3,njitt)*.1 + repmat(grE,1,njitt) ;
                    des = arrayfun(DE,rd(1,:),rd(2,:),rd(3,:));
                    grde = DE(grE(1),grE(2),grE(3));
                    [mn,mni] = min(des);
                    if min(des) > grde
    %                     fprintf('0')
                    else
    %                     fprintf('1')
                        grE = rd(:,mni);
                    end
            end

    %         grEvnp = grE - vn'*diag(vn*grE); %projected into the plan normal to the vertex normal       

            switch constraint
                case 'fn'  % only move along facet planes
                    trv = reshape(dv^-1*repmat(grE,size(dv1,1),1),3,size(dv1,1));        
                    intri = trv(1,:)>0 & trv(2,:)>0; % Find the face that that the gradient leans towards
                    if sum(intri)~=1
                        k
                        warning('Could not identify facet')
                        continue
                    end
                    grEproj = grE' - fn(intri,:)*(fn(intri,:)*grE); 

                case 'vn' % Only move along plane perpendicular to the vertex normal.
                    grEproj = grE' - vn*(vn*grE); 

                case 'none'
                    grEproj = grE'; 
            end

    %           grEprojs = repmat(grE,1,size(dv1,1)) - fn'*diag(fn*grE);       
            Xnew(k,:) = v+grEproj;
            DE = sum((sqrt(sum(dvatt(v+grEproj).^2,2))-d0).^2)-sum((sqrt(sum(dvatt(v).^2,2))-d0).^2);
    %         if DE > 0
    %            keyboard 
    %         end
%             newmesh = TriRep(TRI,Xnew);
              newmesh = TriRep(newmesh.Triangulation,Xnew);
      end


        if mod(ki,10)==0;
            newmeshs{ki/10} = newmesh;
        end
    end
    
         vnl = 0;
         nc = 0;
         while any(vnl < curve_threshold) && nc < 500
            fns = newmesh.faceNormals;
            va = newmesh.vertexAttachments;
            vns=  zeros(0,3);
            for k = 1:length(va), vns(k,:) = mean(fns(va{k},:)); end
            vnl = sqrt(sum(vns.^2,2)); % The length of the vertex normal gives a simple measure of curvature (assuming face normals are pointing out).
            badpts = find(vnl < curve_threshold);
            X = newmesh.X;
%             tri = TRI;
            tri = newmesh.Triangulation;
            nc = nc+1;

            for i = 1:length(badpts)
                 vtr = tri(newmesh.vertexAttachments{badpts(i)},:)';
                 vuq = unique(vtr);
                 vuq(vuq==badpts(i)) = [];
                 X(badpts(i),:) = mean(X(vuq,:));
            end
%             newmesh = TriRep(TRI,X);
            newmesh = TriRep(newmesh.Triangulation,X);
            sum(vnl<curve_threshold)
         end
        
         if any(vnl < curve_threshold)
             %%% Probably means there is a bad connection
             X(badpts,:) = nan;             
             warning('After smoothing %i bad points remained, which will be replaced with nans',sum(badpts));
         end
 end
