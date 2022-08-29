

prp = preop_contact_locations;

if exist('postop_contact_locations')
    pop = postop_contact_locations;
end
%%
if ~exist('pop','var')
    pop=[];
end
mnitr = preop.transforms(3);
rastr  = preop.transforms(1);
if ~exist('blkdat') && exist('lozch')
  blkdat.block='NA';
    blkdat.lozchannels = lozch;
else 
    chlbl = {prp.label};
end

 if exist('blkdat') %|| ~exist('lozch')
    lozch0 = [blkdat.lozchannels,blkdat.hizchannels];
    chlbl = arrayfun(@(x)sprintf('%i: %s %i',x.contact,x.label,x.number),lozch0,'uniformoutput',false);
    lozch = lozch0(ismember(chlbl,{prp.label}));
%     lozch = lozch(ismember(chlbl,regexprep({prp.label},'[.]',':')));
% elseif  ~exist('lozch')
%     error('Need block data')
end
 %%   
fn=sprintf('%s%s%s_contact_locations_fsparc.csv',ddir,filesep,sid);
if exist(fn,'file')
    copyfile(fn,regexprep(fn,'[.]csv','_old.csv'));
end
fid = fopen(fn,'w');

fprintf(fid,'Subject,Channel, Contact, Contact Label,Number,Voxel X,Voxel Y, Voxel Z, RAS X (mm), RAS Y, RAS Z,  MNI X(mm), MNI Y, MNI Z, MNI->mesh X(mm), MNI->mesh Y(mm), MNI->mesh Z(mm), Unwarped X (voxel), Unwarped Y, Unwarped Z');
fprintf(fid,',Desikan-Killiany label,Destrieux label');
fprintf(fid,',,,Channel numbering for block: %s,,,,,Constructed %s by C. Kovach',blkdat.block,datestr(now));

prpi = find(ismember({prp.label},chlbl));
for k = prpi
    
    mnix = mnitr.tr(prp(k).coord);
    rasx = rastr.tr(prp(k).coord);
    mnimeshx = mnitr.tr(t12fs.itr(prxfsadj(k,1:3)));
    %label = sprintf('%s %i',lozch(k).label,lozch(k).number);
%     label = prp(k).label;
    if exist('lozch','var')
        lzch = lozch0(strcmp(chlbl,prp(k).label));
        label = lzch.label;
        fprintf(fid,'\n%s,%i,%i,"%s",%i,',sid,lzch.channel,lzch.contact,label,lzch.number);
    else
        fprintf(fid,'\n%s,%i,%i,"%s",%i,',sid,k,k,prp(k).label,[]);
    end
    fprintf(fid,'%0.2f,',prp(k).coord);
    fprintf(fid,'%0.2f,',rasx);
    fprintf(fid,'%0.2f,',mnix);
    fprintf(fid,'%0.2f,',mnimeshx);
    if ~isempty(pop)
        fprintf(fid,'%0.2f,',pop(k).coord);
    else
        fprintf(fid,'NA,NA,NA,');
    end
    fprintf(fid,'"%s","%s"',dktseg{k},destseg{k}); % Assigned in assign_points_to_Fsparc.m
end
fclose(fid);
    
    

    