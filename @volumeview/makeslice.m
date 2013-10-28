function [SL,edge,trmat,hn] = makeslice(me,nv,pt,rot)

%%% makeslice(me,nv,pt,rot): return a slice normal to nv containing pt rotated about nv by rot


% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 4
    rot = 0;
end


nv = nv./norm(nv);

fig = figure;
set(fig,'visible','off')

sz = size(me.Vol);
sz = sz([2 1 3]);
% h = pcolor((0:sz(2)-1)+2,(0:sz(1)-1)+2,zeros(size(me.Vol(:,:,1)))');
h = pcolor((1:sz(1))+1,(1:sz(2))+1,zeros(size(me.Vol(:,:,1))));
set(h,'zdata',ones(size(me.Vol(:,:,1)))*round(pt(3)+1))

rax = cross( nv,[0 0 1]);
if norm(rax) == 0
    rax = [0 0 1];
    th = 0;
    M = eye(3);
else
    rax = rax./norm(rax);
    th = acos([0 0 1]*nv');
    plax = cross([0 0 1],rax); %axis in plane before rotation
    plax2 = cross(nv,rax);  % axis in plane before rotation
    M = [0 0 1; plax; rax]\[nv;plax2;rax]; % Rotation matrix
end
% axv = [1 0 0; 0 1 0]*M; % axis vectors for the new plane
 rotate(h,rax,-th/pi*180,pt+1);
if rot > 0
    rotate(h,nv,-rot/pi*180,pt+1)
end

 xd = get(h,'xdata');
 yd = get(h,'ydata');
 zd = get(h,'zdata');


hn = slice(me.Vol,yd,xd,zd);

% xy = [xd([1 end]); yd([1 end]); zd([1 end])]';


SL = get(hn,'cdata')';

x0 = [1 1  ;  1 sz(2) ;sz(1) 1 ; sz(1:2) ; 1 1 ];
% x0 = [1 1  ;  1 sz(2) ;sz(1) 1 ; sz(1:2) ];
ln = @(x)x(:);
corn = @(x)ln(x([1 end],[1 end]));

x1 = [corn(xd)-1,corn(yd)-1, corn(zd)-1, ones(4,1)];
x1 = cat(1,x1,x1(1,:) + [nv,0]);
trmat = x1\x0;

if nargout <4
    delete(fig)
else
    set(fig,'visible','on')
end
edge = x1;
