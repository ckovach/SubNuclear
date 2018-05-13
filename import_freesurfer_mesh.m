mnih = readnifti(fullfile(ddir,'..','MNI152_T1_1mm.nii'),true);
% T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))';
fid = fopen(fullfile(ddir,'T1_to_MNI_lin.mat'),'r');
fcont = fread(fid,'uchar=>char')';
fclose(fid);
T = str2num(fcont)';

%%% Directory for freesurfer files for this subject
freesurfer_directory = fullfile(ddir,'..','FS',sprintf('pt%s',sid));

%%%% Get transformation to freesurfer TKR coordinates

[res,out] = system( sprintf('mri_info --vox2ras-tkr %s%smri%sorig.mgz',freesurfer_directory,filesep,filesep));
v2rtk = str2num(out)';
rtk2v = v2rtk^-1;
tktr = transforms('trmat',rtk2v);

%%% Load labels in the lookup table.
[~,out]= system('echo $FREESURFER_HOME');
% fid = fopen('/usr/local/freesurfer/FreeSurferColorLUT.txt','r');
fid = fopen(sprintf('%s/FreeSurferColorLUT.txt',deblank(out)),'r');
txt = fread(fid,'uchar=>char')';
fclose(fid);
txt = regexp(txt,'\n(\d+)\s*([^\s]*)\s*(\d*\s*\d*\s*\d*)[^\n]*','tokens');
txt = cat(1,txt{:});
indx = cellfun(@str2double,txt(:,1));
lut ={};
lut(indx + 1,1) = txt(:,2);
lut(indx + 1,2) = txt(:,3);

%     Heschlloc(i).vox = tktr(i).tr(Heschl
% end
%%


vtkdir = fullfile(freesurfer_directory,'ascii');
d =dir(fullfile(vtkdir,'*.srf'));
if exist(freesurfer_directory,'dir') &&isempty(d)%~exist(vtkdir,'dir')
   warning('No subcortical freesurfer surfaces appear to have been generated for this subject.') 
    if ~exist('aseg2srf','file')
           url='https://dl.dropboxusercontent.com/u/2785709/brainder/2012/aseg2srf';
            fprintf('\nAttempting to download aseg2srf from\n\t%s',url);
            try
                htxt = urlread(fullfile(url,'aseg2srf'));
                fid =fopen('aseg2srf','w');
                fwrite(fid,htxt);
                fclose(fid);
                system('chmod 755 aseg2srf');
            catch
                fprintf('Failed to retrieve aseg2srf.\nCheck\n\t%s\n and try to download manually.\nTerminating...\n',url);
                return
            end
    end
    fprintf('\nRunning aseg2srf. This may take a moment...\n');
    pause(1)
     system(sprintf('./aseg2srf -s "pt%s"',sid))
end
d =dir(fullfile(vtkdir,'*.srf'));

if ~exist(fullfile(ddir,'FST1.nii'),'file')
        res = system( sprintf('mri_convert %s%smri%sorig.mgz %s%sFST1.nii',freesurfer_directory,filesep,filesep,ddir,filesep));
end
fs= volumeview('FST1','FST1.nii',ddir);
clear vvmeshes
for k = 1:length(d)
    
    [~,fn,ext] = fileparts(d(k).name);
    if ~exist(fullfile(vtkdir,[fn,'.vtk']),'file')
        copyfile(fullfile(vtkdir,d(k).name),fullfile(vtkdir,[fn,'.asc']));
        [err,out]=system(sprintf('mris_convert %s.asc %s.vtk',fullfile(vtkdir,fn),fullfile(vtkdir,fn)));
        delete(fullfile(vtkdir,[fn,'.asc']))
    end    
    indx = regexp(fn,'_(\d+)','tokens','once');
    indx = str2double(indx{1})+1;
    lbl = sprintf('FS: %s',lut{indx,1});
    fprintf('Loading %s\n',lbl)
    vtk = readvtk(fullfile(vtkdir,[fn,'.vtk']));

    srf = TriRep(vtk.tri,preop.transforms(1).itr(fs.transforms(1).tr(tktr.tr(vtk.vert))));
    
    vvmesh = preop.addmesh(srf,lbl);
    vvmesh.show=false;
    vvmesh.plotcolor = str2num(lut{indx,2})/255;
    vvmeshes(k)=vvmesh;
end