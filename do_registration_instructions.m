%% Subcortical Parcellation with SubNuclear
%
% 
% 
% The SubNuclear toolbox has been developed with two goals. First is to implement  a relatively straightforward protocol for superimposing atlas-derived substructures on subcortical structures visible in MR's.  The second is to provide a simple and adaptable set of matlab  tools for exploring medical images and localizing electrode positions based imaging artifacts.  
%% Viewing volumes
% 
% 
% SubNuclear is built around the volumeview class.  To instantiate a member of the class, simply type
%%
%  vol = volumeview
%
% at the matlab commandline. You will be prompted to select a volume in .nii or .nii.gz format and the main GUI window will appear with axial, coronal and sagittal views. 
% 
% <</home/ckovach/svn/anatomy_tracing/volumeviewGUI.png>>
% 
% To zoom and pan within the axis windows, use the standard options in the matlab figure menu bar.
%
% To adjust the brightness and contrast, change the intensity_range property of the object. For example,
%
%          vol.intensity_range  = [0 255]; 
%
% adjusts the color scale so that 0 or less is black and 255 or greater is white.
%% Working with volumeview objects
%  
% *1. Passing by reference*
%
% Volumeview objects are handle class objects, which are passed by reference, meaning that if you let
%%
%
%       vol2 = vol
%       vol3 = vol2
%
% vol, vol2 and vol3 now all refer to the same object, so that any change to vol2 will also be reflected in vol3 and vol. If you want to make a copy of an object so that any change of one a leaves the others unaffected, then you must create a new instance of the class with
%
%       vol2 = volumeview(vol)
%
%
% *2. Loading and saving*
%
% 	Volumeview classes can be saved to and loaded from matlab files, however to get the GUI going 
% you need to pass the loaded object as an argument to volumeview:
%
%       load saved_data  vol
%       vol = volumeview(vol)
%
% *3. Methods*
%
% The following functions associated with the volumeview class can be used at the matlab commandline
% 
%       vol.addmesh()  :  Add a mesh to the current volume
%       vol.addpoint()  :  Add a point to the current volume
%       vol.addtransform()  :  add a transform
%       vol.resetaxis()   : reset axes to their original state. Try this if the GUI becomes unresponsive.
%       vol.plotupdate()  :  update the plot.
%
%
% *4. Properties*
%
% The following data can be accessed through the volumeview object
%
%       vol.volumes   :    Medical image class 
%       vol.volume.Data  :   3-d Image data
%       vol.points   :   array of point objects
%       vol.meshes  :  array of mesh objects
%       vol.transforms  :  array of transform objects
%       vol.intensity_range  :  Range used in plotting the volumes.
%       vol.sisters  :  Array of sister volumeview objescts. To make two volumeview objects sisters,  use the command
%
% Sister objects can be linked so that their axis views remain the same. This is especially useful for comparing two volumes.
%
%       vol1.sisters = vol2
% 
% 


sid = inputdlg('Which subject?');
sid= sid{1};

ddir = sprintf('%sfsl.anat',sid);
%% Subcortical parcellation
% 
% For parcellation, SubNuclear uses FSL's module, FIRST, to acquire a gross parcellation of major subcortical structures. These are 
% then coregistered with atlas-derived objects which contain finer-scale structures.  The full sequence of events is 
% implemented in the matlab script, do_registration.  Each section of the script should be run in sequence, though not 
% all at once, i.e. don't run the command do_registration, but evaluate each section in turn. Sections are delimited 
% by' %%'. 
% do_registration first requests the subject ID. All data for a given subject will be placed in the folder 
% [subject_id]fsl.anat
%
%%% FSL
% 
% Subcortical parcellation uses the FIRST module of  FSL (http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIRST/UserGuide), implemented in 
% the script fsl_anat. Obviously this means that FSL has to be installed and in the command line path. FIRST provides a gross 
% parcellation of all the major subcortical structures. Accuracy is not 100%, so in general these parcellations will have to be
% improved with manual editing.  The fsl scripts are self-contained and are generated and run with the matlab command: 
%  

       run_fsl_script

%%%
% You will be prompted to give the preop and postop MR images, which can usually be found on the  HBRL data server  in the 
% .../Images/Image_processing subfolder of a given subject's data folder. Running the complete script takes about an hour on the 
% computational server. While this is running you can begin to choose control points for the warping between preop and postop brains 
% (next section).  
% The result of the script will be:
%
% # The original preop MR in standard orientation ( _sid_fsl.anat/T1.nii.gz)
% # Additional volumes with the  original image linearly and non-linearly  coregistered with the MNI152 standard brain, as well as 
% transforms and displacement fields. 
% # Volumes with labels for parcellated structures.
% # Meshes for each extracted structure, contained in .vtk files.   
% 
%%% Load Data
% 	Meshes can be loaded into the current volume space by right clicking on the mesh list and selecting “Add” or by invoking 
% vol.addmesh at the command line. This will allow you to select a .vtk file, which contains the mesh.  Matlab TriRep objects can also
% be imported using vol.addmesh(TR), where TR is a TriRep object.	After the FSL parcellation is complete, do_registration loads the left 
% and right amygdala and hippocampus  meshes. 
% 
% # Editing Meshes
%    To manually edit the amygdala meshes, first create a copy of the mesh you want to edit then make sure the “Adjust mesh” 
% radio button is enabled. Select the mesh in the mesh list window. Meshes can be edited in three ways:
% # Translation within the plan of a view. To do this click on the mesh in one of the axes and drag it.
% # Expanding the mesh near the current point. Clicking on 'Expand' expands the mesh in the direction normal to the vertex by 
% an amount that depends on the distance to the current point and angle the between the difference vector and the vertex normal. 
% # Similarly 'Contract' has the same effect as 'Expand' only in the opposite direction.

% KTcheck
postop = volumeview(sprintf('%s PostOp',sid),'postop_aligned.nii.gz',ddir);
preop = volumeview(sprintf('%s PreOp',sid),'T1.nii.gz',ddir);
postop.intensity_range = [0 400];
mnih = readnifti('MNI152_T1_1mm.nii',true);

preop.sisters= postop;

%%% Compute the voxel to MNI coordinate transformation
T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))'; %#ok<REMFF1>
tr1 = transforms('trmat',T(1:4,:),'label','T12MNI');
tr2 = transforms('trmat',double(mnih.vox2unit)','label','MNI2mm');
tr = tr1*tr2;
tr.label = 'vox2MNImm';
postop.addtransform(tr);
preop.addtransform(tr);







helpdlg('Check the alignment of the images and select matching control points on the preop and postop brains')
%% Atlas Coregistration for Amygdala
% 
% After you're satisfied with the amygdala parcellation, meshes for amygdalae are loaded from the file,
% mai_template_mni_aligned.mat, and subnuclei are loaded from maiwarp2mni.mat
%  These are extracted from “Atlas of the Human Brain” (Mai 2008). The variables, newmeshR and newmeshL contain
%data for right and left amygdalae, respectively, as TriRep objects  (see matlab help for details on the TriRep class). 
% The vertices in newmeshL and newmeshR are matched to those in the FSL template mesh for the amygdala. Coregistration uses 
% thin plate spline warping to align the vertices in the FSL generated mesh with the atlas-derived mesh. The necessary steps are 
% carried out in do_registration. The resulting warping function is then applied to the vertices of meshes for each subnucleus
% in order to transform the subnuclei into the subject's image space. 
% 

structs = {'L Amyg' , fullfile(ddir,'first_results/T1_first-L_Amyg_first.vtk')
           'R Amyg' , fullfile(ddir,'first_results/T1_first-R_Amyg_first.vtk')};

for i = 1:length(structs)
%     postop.addmesh(structs{i,2},structs{i,1})
     preop.addmesh(structs{i,2},structs{i,1})
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

% preopcontacts = volumeview(preop);
% postopcontacts = volumeview(postop);

% fprintf('Select control points on each image. Make sure they are in the same order')
x1 = cat(1,preop.points.coord); x1(:,4) = 1;
x2 = cat(1,postop.points.coord); x2(:,4) = 1;
% TlinL = x1([1:10,23],:)\x2([1:10,23],:);
% TlinR = x1([7:9,11:22],:)\x2([7:9,11:22],:);

% mx1L = preop.meshes(1).trirep.X; mx1L(:,4) = 1;
% mx2L = postop.meshes(1).trirep.X;mx2L(:,4) = 1;
% TmeshL = mx1L\mx2L;
% % mx2lR = mx1L*TlinR;
% mx1R = preop.meshes(2).trirep.X; mx1R(:,4) = 1;
% mx2R = postop.meshes(2).trirep.X;mx2R(:,4) = 1;
% TmeshR = mx1R\mx2R;
% % mx2lL = mx1R*TlinL;

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
load mai_template_mni_aligned newmeshR newmeshL tmpl2maimatR tmpl2maimatL % Atlas with vertices matcehd to FSL vertices

%%
tr = preop.transforms(3); %This is the vox to MNI transform
% itr = preop.transforms(3).itr(x);


% tr = @(x)x;
% itr = @(x)x;s
warpfunL  = tpswarp(tr.itr(newmeshL.X),preop.meshes(3).trirep.X,.1);  % Left amygdala
% warpfunL =@(x) wrpL(x*diag([-1 1 1]));

warpfunR = tpswarp(tr.itr(newmeshR.X),preop.meshes(4).trirep.X,.1); % Right amygdala

[iwarpfunL,AL]   = tpswarp(tr.tr(preop.meshes(3).trirep.X),newmeshL.X*10);
% iwarpfunL  = @(x) iwpL(x*diag([-1 1 1]));
[iwarpfunR,AR] = tpswarp(tr.tr(preop.meshes(4).trirep.X),newmeshR.X*10);

%%% Affine transformation into template space (so we can match the axial view to the nearest page in the template).
AL(4,4) = 1; AR(4,4) = 1;
tmpl2maimatL(4,4) = 1; tmpl2maimatR(4,4) = 1;
preop.addtransform(tr.trmat*AL*tmpl2maimatL,'Left Amyg. vox 2 atlas');
preop.addtransform(tr.trmat*AR*tmpl2maimatR,'Right Amyg. vox 2 atlas');


kp = find(arrayfun(@(x)~isempty(x.warped),maiwarpL));
kp(end) = [];

colors = {'k','r','g','y','m','c','c','k','b','b'};
for i = 1:length(kp)   
    wrpR = TriRep(maiwarpR(kp(i)).warped.Triangulation,warpfunR(tr.itr(maiwarpR(kp(i)).warped.X)));
    preop.addmesh(wrpR,sprintf('R %s',maiwarpR(kp(i)).label)); 
    
end
for i = 1:length(kp)   
    wrpL = TriRep(maiwarpL(kp(i)).warped.Triangulation,warpfunL(tr.itr(maiwarpL(kp(i)).warped.X)));
    preop.addmesh(wrpL,sprintf('L %s',maiwarpL(kp(i)).label));
end
% 
% for i = 1:length(postopcontacts.points)
% %     contacts(i).coord = postopcontacts.points(i).coord;
% %     contacts(i).label= postopcontacts.points(i).label;
%     
%     postop.addpoint(postopcontacts.points(i).label,imwarp(postopcontacts.points(i).coord));
% end
% 



%% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------