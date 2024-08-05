function [fname,ok]=conn_loadgdfile(fileSharedLink, fileSaveAs, islargefile)
% internal function
% downloads file from google drive

if nargin<3||isempty(islargefile), islargefile=false; end

fileSharedLink=regexprep(fileSharedLink, 'file/d/([^\?\/]*)(.*)$','uc?export=download&id=$1'); % converts file/d/* to uc?export=download&id=* format
[fname,ok]=urlwrite(fileSharedLink,fileSaveAs);

if islargefile 
    try
        hdl=java.io.File(fileSaveAs);
        if hdl.length<1e8
            ok=false;
            str=fileread(fileSaveAs);
            match=regexp(str,'form id\="download-form" action\="([^"]*)"(.*)$','tokens','once');  % confirms GoogleDrive warning message when downloading large files (method 2, NEW)
            if ~ok&&~isempty(match)
                field1=char(regexp(match{2},'name="id"\s*value="([^"]*)"','tokens','once'));
                field2=char(regexp(match{2},'name="uuid"\s*value="([^"]*)"','tokens','once'));
                field3=char(regexp(match{2},'name="confirm"\s*value="([^"]*)"','tokens','once'));
                match=sprintf('%s?id=%s&uuid=%s&export=download&confirm=%s',match{1},field1,field2,field3);
                pause(.25+.5*rand);
                [fname,ok]=urlwrite(match,fileSaveAs);
            end
            match=char(regexp(str,'form id\="download-form" action\="([^"]*)"','tokens','once')); % confirms GoogleDrive warning message when downloading large files (method 1, OLD)
            if ~ok&&~isempty(match)
                match=regexprep(match,'\&amp\;','&');
                pause(.25+.5*rand);
                [fname,ok]=urlwrite(match,fileSaveAs);
            end
        end
    end
end
