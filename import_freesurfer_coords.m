mnih = readnifti(fullfile(ddir,'..','MNI152_T1_1mm.nii.gz'),true);
% T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))';
fid = fopen(fullfile(ddir,'T1_to_MNI_lin.mat'),'r');
fcont = fread(fid,'uchar=>char')';
fclose(fid);
T = str2num(fcont)';

i = 1;
% for i = 1;length(ss)
%     ddir = ss(i).locations.main;
%     T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))'; %#ok<REMFF1> 
%     Tfov = textread(fullfile(ddir,'T1_roi2nonroi.mat'))'; %Transformation into smaller field of view
%     T(1:3,1:3) = abs(diag(diag(preop.transforms(1).trmat(1:3,1:3))))*T(1:3,1:3); % FLIRT uses voxel x size coordinates
%      tr2std = preop.volumes(1).tr2std; % FLIRT also computes transform from standard space;
%     % tr1 = transforms('trmat',Tfov(1:4,:)^-1,'label','tr2fov')*tr2std*transforms('trmat',T(1:4,:),'label','T12MNI');
%     trfov = transforms('trmat',Tfov(1:4,:)^-1,'label','tr2fov'); %Transform to field of view
%     % tr1 = trfov*transforms('trmat',T(1:4,:),'label','T12MNI');
%     tr1 = tr2std*trfov*transforms('trmat',T(1:4,:),'label','T12MNI');
%     tr2 = transforms('trmat',double(mnih.vox2unit)','label','MNI2mm');
%     tr = tr1*tr2;
%     tr.label = 'vox2MNImm';

    freesurfer_directory = fullfile(ddir,'..','FS',sprintf('pt%s',sid));
    [res,out] = system( sprintf('$FREESURFER_HOME/bin/mri_info --vox2ras-tkr %s%smri%sorig.mgz',freesurfer_directory,filesep,filesep));
    v2rtk = str2num(out)';
    rtk2v = v2rtk^-1;
    tktr(i) = transforms('trmat',rtk2v);
%     Heschlloc(i).vox = tktr(i).tr(Heschl
% end
%%

 hemis = 'rl';
% hemis ='l';
fs = [];
clear cxs
prpmesh=meshes;
prpmesh=prpmesh([]);
for k  = 1:length(hemis)
    hemi = hemis(k);
    %   hemi = hemis{k};
    % hemi = 'r';

    cx = read_freesurfer(fullfile(freesurfer_directory,'surf',[hemi,'h.pial']));
    wh = read_freesurfer(fullfile(freesurfer_directory,'surf',[hemi,'h.white']));
%     cx = read_freesurfer(fullfile(freesurfer_directory,'surf',[hemi,'h.orig']));
%     wh = read_freesurfer(fullfile(freesurfer_directory,'surf',[hemi,'h.orig']));

    cxtr = double(tktr(i).tr(cx.vert));
    whtr = double(tktr(i).tr(wh.vert));
    if ~exist(fullfile(ddir,'FST1.nii'))
        res = system( sprintf('$FREESURFER_HOME/bin/mri_convert %s%smri%sorig.mgz %s%sFST1.nii',freesurfer_directory,filesep,filesep,ddir,filesep));
    end
    if isobject(fs) 
      delete(fs)
    end
      fs= volumeview('FST1','FST1.nii',ddir);
    vcx = fs.transforms(1).tr(cxtr);
    vcxs{k} =vcx;
    vwh = fs.transforms(1).tr(whtr);
    vcxtrp = TriRep(cx.tri,double(preop.transforms(1).itr(vcx)));
    %vcxtrp = TriRep(cx.tri,pad(double(preop.transforms(1).itr(vcx)))*Tfov(:,1:3));
    %vcxtrp = TriRep(cx.tri,double(preop.transforms(2).itr(preop.transforms(1).itr(vcx))));
    prpmesh(end+1)=preop.addmesh(vcxtrp,[hemi,' cx']);
    prpmesh(end).show=false;
    prpmesh(end).plotcolor = [0 0 0];
    prpmesh(end).plotargs = {'marker','none','linewidth',1};
    
     vwhtrp = TriRep(wh.tri,double(preop.transforms(1).itr(vwh)));
    % vwhtrp = TriRep(wh.tri,pad(double(preop.transforms(1).itr(vwh)))*Tfov(:,1:3));
    %vwhtrp = TriRep(wh.tri,double(preop.transforms(2).itr(preop.transforms(1).itr(vwh))));
    prpmesh(end+1)=preop.addmesh(vwhtrp,[hemi,' white matter']);
    prpmesh(end).show=false;
    cxs(k).tri = vcxtrp.Triangulation;
    cxs(k).vert = vcxtrp.X;
    
    midgray(k).tri = vcxtrp.Triangulation;
    midgray(k).vert = (vcxtrp.X+vwhtrp.X)/2;
    
    white(k).tri = vwhtrp.Triangulation;
    white(k).vert = vwhtrp.X;
    
end

% 
% if length(cxs)>1
%     cxtrp = TriRep(cat(1,cxs(1).tri,cxs(2).tri+size(cxs(1).vert,1)),preop.transforms(1).tr(cat(1,cxs.vert)));
% else
%     cxtrp = TriRep(cxs(1).tri,preop.transforms(1).tr(cat(1,cxs.vert)));
% end    




