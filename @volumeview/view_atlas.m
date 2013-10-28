% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

function view_atlas
persistent atlas AP tmpl2maimat 


if isempty(atlas)
    load maiatlas2 atlas
    load mai_template_mni_aligned

end


