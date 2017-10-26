function wrpimg= warpVolume(vv,tr,targetvol)

%function warpVolume(vv,tr,sz)
%
% Creates a non-linearly warped volume using the specified transform. 
%
% Inputs: 
%   vv - volumeview object to be transformed
%   tr - function handle for the transformation with IJKtransf = tr(IJK);
%        The transform must be from the target space to the starting space.
%
%   sz - size of destination volume (default is size of starting volume)
%

chunksize = 1e4;

destroy_when_done=false;
if isa(vv,'medimage')
   vv = volumeview(vv);
   destroy_when_done = true;
end
if nargin < 3
    targetvol = vv;
%     sz = size(vv.Vol); % dimensions of input volume.
end

switch class(targetvol)

    case 'medimage'
        wrpimg = medimage(targetvol);
    otherwise
        wrpimg = medimage(targetvol.current.volumes.image);
end



if ~isnumeric(targetvol)
    sz = size(wrpimg.Data);
else
    sz = targetvol;
end


[J,I,K] = meshgrid(1:sz(2),1:sz(1),1:sz(3));

IJK = [I(:),J(:),K(:)]-1;
IJKtr = zeros(size(IJK));
%%
% tic
h = waitbar(0,'Applying transform...');
for i = 1:chunksize:size(IJK,1)
%     fprintf('\nDoing points %i through %i of %i...',i,chunksize+i-1,size(IJK,1));
    IJKtr(i:min(chunksize-1+i,end),:)=tr(IJK(i:min(chunksize-1+i,end),:)) ;
    waitbar(i/size(IJK,1),h)
%     i
%    toc 
end
delete(h)

vol = reshape(vv.current.volumes.interpolant(IJKtr),sz);

wrpimg.Data = vol;
wrpimg.file = ['warped_',wrpimg.file];

if destroy_when_done
    delete(vv);
end
