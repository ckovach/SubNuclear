
function Tr = maketransform(vfr,vto)

flds = {'meshes','points'};
inp = listdlg('liststring',flds,'promptstring','Make transform from which data?','selectionmode','single');

fld = flds{inp};

go = 1;
while go
    inpfr = listdlg('liststring',{vfr.(fld).label},'promptstring','Choose "from" points.','selectionmode','multi');
    inpto = listdlg('liststring',{vto.(fld).label},'promptstring','Choose "from" points.','selectionmode','multi');
    if isempty(inpfr) || isempty(inpto)
        return
    end
    switch fld
        case 'meshes'
            Xfr = cat(1,vfr.(fld)(inpfr).trirep.X);
            Xto = cat(1,vto.(fld)(inpto).trirep.X);
        case 'points'
            Xfr = cat(1,vfr.(fld)(inpfr).coord);
            Xto = cat(1,vto.(fld)(inpto).coord);        
    end

    go = size(Xfr,1)~=size(Xto,1);
    if go
        warning('Number of points doesn''t match')
        h =warndlg('Number of points doesn''t match. Try again.','','modal');
        uiwait(h);
        
    end
end

Xfr(:,4) = 1;
Xto(:,4) = 1;

T = Xfr\Xto;

if nargout >0
    Tr = transforms('trmat',T);
else

    vto.addtransform(T)
    if vfr ~=vto
        vfr.addtransform(vto.transforms(end));
    end
end