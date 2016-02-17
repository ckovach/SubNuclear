

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%% Get Subject ID
sid = inputdlg('Which subject?');
sid= sid{1};

ddir = sprintf('%sfsl.anat',sid);

%% Subcortical parcellation
% Calls FSL to 
run_fsl_script

%% Load Data

% KTcheck('volumeview')
preop = volumeview(sprintf('%s PreOp',sid),'T1_fullfov.nii.gz',ddir);
postop = volumeview(sprintf('%s PostOp',sid),'postop_aligned.nii.gz',ddir);
postop.intensity_range = [0 800];

preop.sisters= postop;
helpdlg('Make sure that preop and postop images are correctly aligned. If they aren''t than realign with FLIRT.')

%% Compute the voxel to MNI coordinate transformation

%%% The parameters in the fsl mat transform files are idiosyncratic so best
%%% to compute the transform with fsl's img2imgcoord utility.

[res,out] = system('echo $FSLDIR');
if ~isempty(out)
    fsl_imagedir = fullfile(deblank(out),'..','data','standard');
    if ~exist(fsl_imagedir,'dir')
        fsl_imagedir = fullfile(deblank(out),'data','standard');
    end
    fsl_imagedir = regexprep(fsl_imagedir,'[\n]','');
else
    fsl_imagedir = '.';
end

tr2std = preop.volumes(1).tr2std; % FLIRT also computes transform from standard space;

com = 'printf ''0 0 0\n1 0 0\n0 1 0\n 0 0 1'' | img2imgcoord -vox -src %s -dest %s -xfm %s';
A = [0 0 0 1; 1 0 0 1; 0 1 0 1; 0 0 1 1]; 
[res,out] = system(sprintf(com,fullfile(ddir,'T1.nii.gz'),...
                         fullfile(ddir,'T1_fullfov.nii.gz'),...
                         fullfile(ddir,'T1_roi2nonroi.mat')));
out = regexprep(out,'.*[:)]','');
xfov = str2num(out);
xfov(:,end+1) = 1;
Tfov = xfov\A; %%% T1_fullfov to T1 transform

%%% Postop to preop linear transform
% [res,out] = system(sprintf(com,fullfile(ddir,'postop_orig.nii.gz'),...
%                          fullfile(ddir,'T1_fullfov.nii.gz'),...
%                          fullfile(ddir,'post_to_pre.mat')));
%  if res==0                    
%     out = regexprep(out,'.*:','');
%     poapt = str2num(out);
%     poapt(:,end+1) = 1;
%     Tpo2poa = poapt\A; %%% Postop orig to preop linear
%  end
 
[res,out] = system(sprintf(com,fullfile(ddir,'T1.nii.gz'),...
                         fullfile(fsl_imagedir,'MNI152_T1_1mm.nii.gz'),...
                         fullfile(ddir,'T1_to_MNI_lin.mat')));
out = regexprep(out,'.*[:)]','');

xmni = str2num(out);
xmni(:,end+1) = 1;

mnih = readnifti(fullfile(fsl_imagedir,'MNI152_T1_1mm.nii.gz'),true);
Tmni =Tfov*(A\xmni)*mnih.vox2unit';
tr = transforms('trmat',Tmni,'label','vox2MNImm');
% mnih = readnifti('MNI152_T1_1mm.nii',true);
% T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))'; %#ok<REMFF1> 
% Tfov = textread(fullfile(ddir,'T1_roi2nonroi.mat'))'; %Transformation into smaller field of view
% T(1:3,1:3) = diag(max(abs(preop.transforms(1).trmat(1:3,1:3))))*T(1:3,1:3); % FLIRT uses voxel x size coordinates
% Tfov(4,1:3) =Tfov(4,1:3)* diag(max(abs(preop.transforms(1).trmat(1:3,1:3))).^-1);
% tr2std = preop.volumes(1).tr2std; % FLIRT also computes transform from standard space;
% % tr1 = transforms('trmat',Tfov(1:4,:)^-1,'label','tr2fov')*tr2std*transforms('trmat',T(1:4,:),'label','T12MNI');
% trfov = transforms('trmat',Tfov(1:4,:)^-1,'label','tr2fov'); %Transform to field of view
% % tr1 = trfov*transforms('trmat',T(1:4,:),'label','T12MNI');
% tr1 = tr2std*trfov*transforms('trmat',T(1:4,1:3),'label','T12MNI');
% tr2 = transforms('trmat',double(mnih.vox2unit)','label','MNI2mm');
% tr = tr1*tr2;
tr.label = 'vox2MNImm';
postop.addtransform(tr);
preop.addtransform(tr);

[res,out] = system(sprintf(com,fullfile(ddir,'postop_orig'),...
                         fullfile(ddir,'T1_fullfov'),...
                         fullfile(ddir,'post_to_pre.mat')));
out = regexprep(out,'.*[:)]','');
xp2p= str2num(out);
xp2p(:,end+1) = 1;
Tpo2pr  = A\xp2p;

%save(fullfile(ddir,sprintf('%s_volview',sid)),'preop','postop'),



helpdlg('Check the alignment of the images and select matching control points on the preop and postop brains')
%% Atlas Coregistration for Amygdala
% 
% After you're satisfied with the amygdala parcellation, meshes for amygdalae are loaded from the 
%file,  mai_template_mni_aligned.mat, and subnuclei are loaded from maiwarp2mni.mat
%  These are extracted from “Atlas of the Human Brain” (Mai 2008). The variables, newmeshR and 
%newmeshL contain data for right and left amygdalae, respectively, as TriRep objects  (see matlab 
%help for details on the TriRep class). The vertices in newmeshL and newmeshR are matched to those
%in the FSL template mesh for the amygdala. Coregistration uses thin plate spline warping to align 
%the vertices in the FSL generated mesh with the atlas-derived mesh. The necessary steps are carried
%out in do_registration. The resulting warping function is then applied to the vertices of meshes
%for each subnucleus in order to transform the subnuclei into the subject's image space. 
% 

structs = {'L Amyg' , fullfile(ddir,'first_results/T1_first-L_Amyg_first.vtk')
           'R Amyg' , fullfile(ddir,'first_results/T1_first-R_Amyg_first.vtk')
           'L Hipp' , fullfile(ddir,'first_results/T1_first-L_Hipp_first.vtk')
           'R Hipp' , fullfile(ddir,'first_results/T1_first-R_Hipp_first.vtk')};
%%% The vtk file contain coordinates for the brain after it has been
%%% reoriented into RAS space and put into a smaller FOV.
%%% This still needs to be verified with non-standard orientations.
trfov = transforms('trmat',Tfov(1:4,1:3),'label','tr2fov'); %Transform to field of view
preop2fsl = preop.volumes.tr2std...    %    *trfov...
    *transforms('trmat',Tfov(1:4,:)*diag([max(abs(preop.transforms(1).trmat(1:3,1:3))) 1])*eye(4,3));
       
for i = 1:length(structs)
%     postop.addmesh(structs{i,2},structs{i,1})
     vtk(i) = readvtk(structs{i,2});
     trp = TriRep(vtk(i).tri,preop2fsl.itr(vtk(i).vert));
     preop.addmesh(trp,structs{i,1})
end

helpdlg('Now check and if necessary adjust the amygdalae meshes')


%%
% helpdlg('Now adjust the meshes')
% 
% for i = 3:length(preop.meshes)
%    postop.addmesh(preop.meshes(i));
%    
% end


%% Save progress
save(fullfile(ddir,sprintf('%s_volview',sid)),'preop','postop'),
%% Next apply TPS warping


hd = helpdlg('Select control points on each image. Make sure they are in the same order');
uiwait(hd)
x1 = cat(1,preop.current.points.coord); x1(:,4) = 1;
x2 = cat(1,postop.current.points.coord); x2(:,4) = 1;



mwarp = tpswarp(x1(:,1:3),x2(:,1:3)); % this warps from preop to postop
imwarp = tpswarp(x2(:,1:3),x1(:,1:3)); % this from postop to preop
% % mx2tps  = mwarp(mx1(:,1:3));
% 
% XL = preop.meshes(end-1).trirep.X; XL(:,4)=1;
% % mlinL = TriRep(preop.meshes(1).trirep.Triangulation,XL*TlinL(:,1:3));
% mtpsL = TriRep(preop.meshes(end-1).trirep.Triangulation,mwarp(XL(:,1:3)));
% XR = preop.meshes(end).trirep.X; XR(:,4)=1;
% % mlinR = TriRep(preop.meshes(2).trirep.Triangulation,XR*TlinR(:,1:3));
% mtpsR = TriRep(preop.meshes(end).trirep.Triangulation,mwarp(XR(:,1:3)));


%%
figure(postop.fig)
helpdlg('Now add the contact locations')
%% Project contacts to preop brain
figure(postop.fig)
hd = helpdlg('Select contact points in the postop brain, and press OK');
uiwait(hd)
for i =1:length(postop.current.points);
    pt= postop.current.points(i);
    
    preop.addpoint(pt.label,imwarp(pt.coord));
end

% %%  Project the mesh from the preop brain onto the postop as an additional sanity check
% % postop.addmesh(mlinL,'LA Lin');
% postop.addmesh(mtpsL,'LA TPS');
% % postop.addmesh(mlinR,'RA Lin');
% postop.addmesh(mtpsR,'RA TPS');
% 
%% Now load the mai atlas subnuclei
load maiwarp2mni  % Atlas warped into MNI space based on FSL meshes
load mai_template_mni_aligned newmeshR newmeshL% tmpl2maimatR tmpl2maimatL % Atlas with vertices matcehd to FSL vertices

%%
tr = preop.transforms(end); %This sould be the vox to MNI transform

figure(preop.fig)
h = msgbox('Select the Left Amygdala mesh and press OK');
uiwait(h)
drawnow
warpfunL  = tpswarp(tr.itr(newmeshL.X),preop.current.meshes.trirep.X,.01);  % Left amygdala
[iwarpfunL,AL]   = tpswarp(tr.tr(preop.current.meshes.trirep.X),newmeshL.X,.01);

figure(preop.fig)
h = msgbox('Select the Right Amygdala mesh and press OK');
uiwait(h)
drawnow
warpfunR = tpswarp(tr.itr(newmeshR.X),preop.current.meshes.trirep.X,.01); % Right amygdala
[iwarpfunR,AR] = tpswarp(tr.tr(preop.current.meshes.trirep.X),newmeshR.X,.01);

%%% Affine transformation into template space (so we can match the axial
%%% view to the nearest page in the template).
AL(4,4) = 1; AR(4,4) = 1;
tmpl2maimatL(4,4) = 1; tmpl2maimatR(4,4) = 1;
preop.addtransform(tr.trmat*AL*tmpl2maimatL,'Left Amyg. vox 2 atlas');
preop.addtransform(tr.trmat*AR*tmpl2maimatR,'Right Amyg. vox 2 atlas');


kp = find(arrayfun(@(x)~isempty(x.warped),maiwarpL));
kp(end) = [];

colors = {'k','r','g','y','m','c','c','k','b','b'};

lineprops = {'k','r','g','y','m','k','k','k','b',[1 .5 .5]};
for i = 1:length(kp)   
    wrpR = TriRep(maiwarpR(kp(i)).warped.Triangulation,warpfunR(tr.itr(maiwarpR(kp(i)).warped.X)));
    preop.addmesh(wrpR,sprintf('R %s',maiwarpR(kp(i)).label)); 
    preop.meshes(end).plotcolor = lineprops{i};
    preop.meshes(end).plotargs = {'linewidth',2};
end
for i = 1:length(kp)   
    wrpL = TriRep(maiwarpL(kp(i)).warped.Triangulation,warpfunL(tr.itr(maiwarpL(kp(i)).warped.X)));
    preop.addmesh(wrpL,sprintf('L %s',maiwarpL(kp(i)).label));
    preop.meshes(end).plotcolor = lineprops{i};
    preop.meshes(end).plotargs = {'linewidth',2};
end

save(fullfile(ddir,sprintf('%s_volview',sid)),'preop','postop'),


%% Conversion to Freesurfer tkRAS
%
% This checks whether freesurfer bin files are in the OS path, then adds
% transforms to freesurfer space to the volumeview if they are.
%

freesurfer_path = fullfile(pwd,'FS');  % Path to freesurfer subject folder  

%%% Uncomment for defailt FS data path
%[res,out] = system('echo $FREESURFER_HOME');
%freesurfer_path = fullfile(deblank(out),'data');

freesurfer_files = fullfile(freesurfer_path,sprintf('pt%s',sid));

[res,out] = system( sprintf('mri_convert %s%smri%sorig.mgz %s%sFST1.nii',freesurfer_files,filesep,filesep,ddir,filesep));
if res~=0
    warning(['Error encountered while attempting to import the subject %s''s volume from freesurfer:\n\n\t%s\n\nIf you want to interface with freesurfer ',...
              'make sure the appropriate binaries\nare functioning and in the current system OS path, and that required data',...
              '\nare located at\n\t%s\n\n...otherwise,ignore this warning.\n\n'],sid,out,freesurfer_directory);
else
    fs = volumeview('FST1','FST1.nii',ddir);
    [res,out] = system( sprintf('mri_info --vox2ras-tkr %s%smri%sorig.mgz',freesurfer_files,filesep,filesep));
     v2rtk = str2num(out)';
    v2fs = preop.transforms(1).trmat*fs(1).transforms(1).trmat^-1;
    v2fsras = v2fs*v2rtk;

    preop.addtransform(v2fs(:,1:3),'freesurfer vox');
    preop.addtransform(v2fsras(:,1:3),'freesurfer tkRAS');
    postop.addtransform(v2fs(:,1:3),'freesurfer vox');
    postop.addtransform(v2fsras(:,1:3),'freesurfer tkRAS');

end

%% Import freesurfer meshes
