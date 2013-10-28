function A = ReadTabDelim(fname,delim)

% A revolutionary new program that cleverly takes the tab delimited contents of
% a file and reads it into a cell in string format. 

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

fid = fopen(fname);
A = {};
i = 1;

if nargin < 2
    delim = '\t';
end

% fprintf('\n\nline       ');
while ~feof(fid)
    
    s = fgetl(fid);
%     j = 1;
    
    a = regexp(s,sprintf('([^%s]*)%s?',delim,delim),'tokens');
    if ~isempty(a)
        A(i,1:length(a)) = [a{:}];
    end

    
    i = i+1;
%     fprintf('\b\b\b\b\b\b%06i',i) ;
end

fclose(fid);