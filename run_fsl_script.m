
%Generates and, optionally, runs FSL scripts to coregister and parcellate
%MR images.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

run_script = false; %generates and runs script if true, otherwise generates script only
% 
% srvdir = '/hbrl/Data/';
% if ~exist(srvdir,'dir')
srvdir = pwd;
% end
odir = sprintf('%sfsl',sid);
if exist([odir,'.anat'],'dir')
     yn = questdlg(sprintf('The directory %s already exists. Running the FSL script will erase ALL contents. Rename the original directory if you do not want that to happen.\nContinue?',[odir,'.anat']));
     if ~isequal(yn,'Yes')
         return
     end
end

[fnpreop,pthpreop] = uigetfile(fullfile(srvdir,'*.nii*'),'Load PREOP image');
[fnpostop,pthpostop] = uigetfile(fullfile(pthpreop,'*.nii*'),'Load POSTOP image');
 
scrfname = sprintf('run_fsl_%s.sh',sid); % bash script name

[pth,fn,ext] = fileparts(fnpreop);
if isequal(ext,'.gz')
    ext = '.nii.gz';
end
com1 = sprintf('cp %s %sT1temp_orig%s',fullfile(pthpreop,fnpreop),sid,ext);
com2 = sprintf('fslreorient2std %sT1temp_orig  %sT1temp',sid,sid);
com3 = 'echo ''Linear coregistration of preop and postop images...''';
com4 = sprintf('flirt -in %s -ref %sT1temp -o %spostop_aligned -omat %spost_to_pre.mat  -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp spline &',fullfile(pthpostop,fnpostop),sid,sid,sid);
com5 = sprintf('fsl_anat --weakbias  -i %sT1temp -o %s --noreorient --clobber',sid,odir);
com6 =  sprintf('fslmaths %spostop_aligned  %s.anat%spostop_aligned',sid,odir,filesep);
com7 =  sprintf('fslmaths %s %s.anat%spostop_orig',fullfile(pthpostop,fnpostop),odir,filesep);
com8 =  sprintf('mv %spost_to_pre.mat  %s.anat%spost_to_pre.mat',sid,odir,filesep);
com9 = sprintf('rm %sT1temp* %spostop_aligned.mat',sid,sid);
% com5 = sprintf('fsl_anat --weakbias  -d %sfsl -o %sfsl --noreorient',sid,sid);
% com2 = sprintf('flirt -in %s -ref %sfsl.anat/T1.nii.gz -o %sfsl.anat/%s -omat %s/post_to_pre  -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear',fullfile(pthpostop,fnpostop),sid,sid,'post_op_aligned',odir);
% com1 = sprintf('fsl_anat --weakbias  -i %s -o %sfsl -noreorient',fullfile(pthpreop,fnpreop),sid);





fid = fopen(scrfname,'w');

fprintf(fid,'#!/bin/bash\n\n');
fprintf(fid,'printf ''Running fsl commands:\\n\\n''\n');
fprintf(fid,'cat %s\n',scrfname);
fprintf(fid,'printf ''\\n\\n''\n\n');
% fprintf(fid,com1);
fprintf(fid, '\n\n%s',com1,com2, com3,com4,com5,com6,com7,com8,com9);
% fprintf(fid,com2);

fclose(fid);
system(sprintf('chmod 755 %s',scrfname)) ;
%% Now run the script in a separate terminal
if run_script
    [stat,msg] = system(sprintf('xterm -e "./%s;read" &',scrfname) );
else
    fprintf('\nRun the FSL pipeline with the command:\n\n\t[stat,msg] = system(''xterm -e "./%s;read" &'')\n\n',scrfname)
end