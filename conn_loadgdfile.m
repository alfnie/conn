function [fname,ok]=conn_loadgdfile(fileSharedLink, fileSaveAs, islargefile)
% internal function
% downloads file from google drive

if nargin<3||isempty(islargefile), islargefile=false; end

fileSharedLink=regexprep(fileSharedLink, 'file/d/([^\?\/]*)(.*)$','uc?export=download&id=$1'); % converts file/d/* to uc?export=download&id=* format
[fname,ok]=urlwrite(fileSharedLink,fileSaveAs);

if islargefile 
    try
        hdl=java.io.File(fileSaveAs);
        if hdl.length<1e9
            str=fileread(fileSaveAs);
            match=char(regexp(str,'form id\="download-form" action\="([^"]*)"','tokens','once')); % confirms GoogleDrive warning message when downloading large files
            if ~isempty(match)
                match=regexprep(match,'\&amp\;','&');
                [fname,ok]=urlwrite(match,fileSaveAs);
            end
        end
    end
end
