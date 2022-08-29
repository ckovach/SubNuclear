function [err,res] = pdflink(fig,txt,hyperlinks,fname,varargin)

%   pdflink(fig,txt,hyperlinks,fname)
%
% This script prints a figure to pdf with hyperlinks attached to the
% specified text objects. 
%
% INPUTS:
%
%   fig  - figure handle
%   txt  - array of text objects contained within the figure.
%   hyperlinks - cell array of hyperlinks assigned to corresponding txt
%               array objects.
%   fnanme  -  name of the output file. 
%
%
% This script make external calls to ghostwrite, which must be available at
% the command line for this to work.


% C, Kovach 2018

if ischar(hyperlinks)
    hyperlinks = {hyperlinks};
end

if length(hyperlinks)~=length(txt)
   error('Number of hyperlinks needs to match the number of text objects') 
end

switch computer
    case {'GLNXA32','GLNXA64','MACI32','MACI64'}
        gs = 'gs';
    case 'PCWIN64'
        gs='"C:\Program Files\gs\gs9.22\bin\gswin64c"';
    case 'PCWIN32'
        gs='"C:\Program Files\gs\gs9.22\bin\gswin32c"';
end
    
res = system([gs,' --version']);

if res > 0
    gslink = 'https://www.ghostscript.com/download/gsdnld.html';
    error('%s calls ghostscript.\nPlease make sure ghostscript is installed and available at the command line.\nFor details, visit \n\t\t%s',upper(mfilename),gslink)
end
unit = 'points';
set(fig,'Units',unit,'PaperUnits',unit);
set(txt,'Units',unit);

ppos = get(fig,'paperposition');
set(fig,'papersize',ppos(3:4),'paperpos',[0 0 ppos(3:4)],'paperpositionmode','auto');



gscom =...
[' -c "[/Rect [%i %i %i %i]"',...
' -c "  /Border [0 0 1]"',...   %Change hyperlink bounding box thickness here
' -c "  /Color [.7 0 0]"',...  %Change hyperlink bounding box color here
' -c "  /Page 1"',...
' -c "  /Action <</Subtype /URI"',...
' -c "  /URI (%s)>>" ',...
' -c "  /Subtype /Link"',... 
' -c "  /ANN pdfmark"'];


gscoms = '';         
gscomfile = tempname;
fid = fopen(gscomfile,'w');

for k = 1:length(txt)
    
    ax = get(txt(k),'parent');
    set(ax,'Units','points');

    gt = txt(k);
    ext = get(gt,'Extent');
    xpos = get(gt,'Position');
    ppos = get(fig,'PaperPosition');
    pos = ext([1 2 1 2])+[0 0 ext([ 3 4])] +xpos([1 2 1 2])+ppos([1 2 1 2]);
    gscoms = fprintf(fid,['\n',gscom],round(pos),hyperlinks{k});
%     gscoms = sprintf(['%s',gscom],gscoms,round(pos),hyperlinks{k});
end
fclose(fid);

[pth,fn]=fileparts(fname);
tmp = tempname;

print(fig,varargin{:},'-dpsc',tmp);
[err,res]=system(sprintf('%s -sDEVICE=pdfwrite -dEPScrop -o %s.pdf @%s -f %s.ps ',...
                               gs,fullfile(pth,fn),gscomfile,tmp));
                           
