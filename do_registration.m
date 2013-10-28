

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
postop.intensity_range = [0 400];

preop.sisters= postop;
helpdlg('Make sure that preop and postop images are correctly aligned. If they aren''t than realign with FLIRT.')

%% Compute the voxel to MNI coordinate transformation
mnih = readnifti('MNI152_T1_1mm.nii',true);
T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))'; %#ok<REMFF1> 
Tfov = textread(fullfile(ddir,'T1_roi2nonroi.mat'))'; %Transformation into smaller field of view
T(1:3,1:3) = abs(diag(diag(preop.transforms(1).trmat(1:3,1:3))))*T(1:3,1:3); % FLIRT uses voxel x size coordinates
 tr2std = preop.volumes(1).tr2std; % FLIRT also computes transform from standard space;
% tr1 = transforms('trmat',Tfov(1:4,:)^-1,'label','tr2fov')*tr2std*transforms('trmat',T(1:4,:),'label','T12MNI');
trfov = transforms('trmat',Tfov(1:4,:)^-1,'label','tr2fov'); %Transform to field of view
% tr1 = trfov*transforms('trmat',T(1:4,:),'label','T12MNI');
tr1 = tr2std*trfov*transforms('trmat',T(1:4,:),'label','T12MNI');
tr2 = transforms('trmat',double(mnih.vox2unit)','label','MNI2mm');
tr = tr1*tr2;
tr.label = 'vox2MNImm';
postop.addtransform(tr);
preop.addtransform(tr);







helpdlg('Check the alignment of the images and select matching control points on the preop and postop brains')
%% Atlas Coregistration for Amygdala
% 
% After you're satisfied with the amygdala parcellation, meshes for amygdalae are loaded from the file,  mai_template_mni_aligned.mat, and subnuclei are loaded from maiwarp2mni.mat
%  These are extracted from “Atlas of the Human Brain” (Mai 2008). The variables, newmeshR and newmeshL contain data for right and left amygdalae, respectively, as TriRep objects  (see matlab help for details on the TriRep class). The vertices in newmeshL and newmeshR are matched to those in the FSL template mesh for the amygdala. Coregistration uses thin plate spline warping to align the vertices in the FSL generated mesh with the atlas-derived mesh. The necessary steps are carried out in do_registration. The resulting warping function is then applied to the vertices of meshes for each subnucleus in order to transform the subnuclei into the subject's image space. 
% 

structs = {'L Amyg' , fullfile(ddir,'first_results/T1_first-L_Amyg_first.vtk')
           'R Amyg' , fullfile(ddir,'first_results/T1_first-R_Amyg_first.vtk')
           'L Hipp' , fullfile(ddir,'first_results/T1_first-L_Hipp_first.vtk')
           'R Hipp' , fullfile(ddir,'first_results/T1_first-R_Hipp_first.vtk')};
%%% The vtk file contain coordinates for the brain after it has been
%%% reoriented into RAS space and put into a smaller FOV.
preop2fsl = preop.volumes.tr2std*trfov;
       
for i = 1:length(structs)
%     postop.addmesh(structs{i,2},structs{i,1})
     vtk(i) = readvtk(structs{i,2});
     trp = TriRep(vtk(i).tri,preop2fsl.itr(vtk(i).vert));
     preop.addmesh(trp,structs{i,1})
end

helpdlg('Now adjust the amygdalae meshes')


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


helpdlg('Select control points on each image. Make sure they are in the same order')
x1 = cat(1,preop.current.points.coord); x1(:,4) = 1;
x2 = cat(1,postop.current.points.coord); x2(:,4) = 1;



mwarp = tpswarp(x1(:,1:3),x2(:,1:3)); % this warps from preop to postop
imwarp = tpswarp(x2(:,1:3),x1(:,1:3)); % this from postop to preop
% mx2tps  = mwarp(mx1(:,1:3));

XL = preop.meshes(end-1).trirep.X; XL(:,4)=1;
% mlinL = TriRep(preop.meshes(1).trirep.Triangulation,XL*TlinL(:,1:3));
mtpsL = TriRep(preop.meshes(end-2).trirep.Triangulation,mwarp(XL(:,1:3)));
XR = preop.meshes(end).trirep.X; XR(:,4)=1;
% mlinR = TriRep(preop.meshes(2).trirep.Triangulation,XR*TlinR(:,1:3));
mtpsR = TriRep(preop.meshes(end-1).trirep.Triangulation,mwarp(XR(:,1:3)));


%%
figure(postop.fig)
helpdlg('Now add the contact locations')
%% Project contacts to preop brain
for i =length(preop.points)+1:length(postop.points);
    pt= postop.points(i);
    
    preop.addpoint(pt.label,imwarp(pt.coord));
end

%%  Project the mesh from the preop brain onto the postop as an additional sanity check
% postop.addmesh(mlinL,'LA Lin');
postop.addmesh(mtpsL,'LA TPS');
% postop.addmesh(mlinR,'RA Lin');
postop.addmesh(mtpsR,'RA TPS');

%% Now load the mai atlas subnuclei
load maiwarp2mni  % Atlas warped into MNI space based on FSL meshes
load mai_template_mni_aligned newmeshR newmeshL% tmpl2maimatR tmpl2maimatL % Atlas with vertices matcehd to FSL vertices

%%
tr = preop.transforms(end); %This sould be the vox to MNI transform

figure(preop.fig)
h = msgbox('Select the Left Amygdala mesh and press OK');
uiwait(h)
warpfunL  = tpswarp(tr.itr(newmeshL.X),preop.current.meshes.trirep.X,.01);  % Left amygdala
[iwarpfunL,AL]   = tpswarp(tr.tr(preop.current.meshes.trirep.X),newmeshL.X,.01);

figure(preop.fig)
h = msgbox('Select the Right Amygdala mesh and press OK');
uiwait(h)
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

