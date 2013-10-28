function kt_present = KTcheck(check_file)

% This function checks if KovachTools (KT) is present in the current path. If
% not, it offers to attempt a fresh installation through svn. It returns
% a value of 1 if GR is present or successfully installed and 0 otherwise.
%
% Subversion can be downloaded at http://subversion.apache.org/packages.html
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if nargin < 1
    check_file = 'svn_add_all_files.m';
else
    [pth,fn] = fileparts(check_file); %#ok<*ASGLU>
    check_file =[fn,'.m'];
end


%%% URL for current KovachTools repository
grurl = 'https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/';

if ~exist(check_file,'file')
    beep
    
    
   fprintf(['\n--------MISSING FILES---------\nNecessary files appear not to be in the current path.\n\nPlease install KovachTools and add to path\n',...
            '%s\n\n',grurl]) 
    inp = 'x';    
    while ~ismember(lower(inp),{'y','n'})
        inp = input(sprintf(['\nDo you want me to try to install KovachTools now?\n(This requires',...
                    ' command line svn, and username/passwd access)\nY/N:']),'s');
    end
    
    if isequal(lower(inp),'y')
        
        
        fprintf('\nChoose a place to install..')
        installdir = '';
        while exist(fullfile(installdir,'.svn'),'dir') > 0
           installdir = fullfile('..',installdir); 
        end
                
        installdir = uigetdir(installdir,'Choose where to install KovachTools.');
        while exist(fullfile(installdir,'.svn'),'dir') > 0  && ~exist(fullfile(installdir,'misctools'),'dir')
            warndlg('Selected directory must not be a working copy of a subversion repository. Please choose another.')
            installdir = uigetdir(installdir,'Selected directory must not be a subversion repository.');        
            if installdir == 0
                return
            end
        end    
        
        if exist(installdir,'dir') > 0 && exist(fullfile(installdir,'misctools'),'dir')
            fprintf('\nAdding existing version to the matlab path.')
            addpath(installdir)
            kt_present = true;
            return
        elseif exist(installdir,'dir') > 0
            grpath = fullfile(installdir,'KovachTools');
        else
            grpath = installdir;
        end
        
        com = sprintf('svn co %s %s',grurl,grpath);
        
        %%% Attempt to checkout using subversion at the system command-line
        [stat,res] = system(com);  
        
        if stat >0  %%% If command returned an error
            beep
            if ispc                  
                    url = sprintf('download a free version at\n\n\t\thttp://www.sliksvn.com/en/download');
            elseif ismac 
                    url = sprintf('download a free version at\n\n\t\thttp://www.open.collab.net/downloads/community/');
            elseif isunix
                    url = sprintf('install it from your \nrepository (eg run  ''sudo apt-get install svn'')');
            else
                    url = sprintf('find\na version suitable for your system at\n\n\thttp://subversion.apache.org/packages.html');                    
            end
            fprintf(['\n Installation failed with error \n\n\t%s.\n\nIf this happened because you don''t have a working',...
                     ' copy of SVN, you can %s\n\n'],res,url)

        else
            addpath(grpath)
        end
        
    else
        stat = 1;
    end
    
    
else
    pth = fileparts(which(check_file));
    fprintf('\nUsing KovachToolbox at: \n\n')
    [statwc,res]=system(sprintf('svn info %s',pth));
    fprintf(res);
    if statwc
        warning('A necessary file was found, but not within a working copy of KovachToolbox.')
    end        
    stat = 0;
    
end

kt_present = stat==0;