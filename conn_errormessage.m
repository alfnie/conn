function [str,PrimaryMessage]=conn_errormessage(me,projectname,dodebug,connver,rwhos)

if nargin<2, projectname=[]; end
if nargin<3||isempty(dodebug), dodebug=0; end
if nargin<4||isempty(connver), connver=conn('ver'); end
if nargin<5, rwhos=[]; end

stradd={};
PrimaryMessage='display full error message';
if str2num(datestr(version('-date'),'yyyy.mm'))>=2008.08 %Matlab>=2008b
    str=regexprep(char(getReport(me,'basic','hyperlinks','off')),'[\t ]+',' ');
    if strfind(str,'<DisregardMessage>'), str='Process terminated by user'; PrimaryMessage=str;
    else
        idx0=strfind(str,'<nodetailsflag>');
        if ~isempty(idx0),
            str=str(idx0(1)+numel('<nodetailsflag>'):end);
        else
            str=regexprep(char(getReport(me,'extended','hyperlinks','off')),'[\t ]+',' ');
            if dodebug
                try
                    str1=rwhos;
                    for n=1:numel(str1),
                        try
                            if ~ismember(str1(n).name,{'PrimaryMessage','str','stradd','ans','boffset','me'})
                                str2=evalc(['disp(' str1(n).name ')']);
                                str2=str2(str2>=32);
                                if numel(str2)>200, str2=[str2(1:min(200,numel(str2))),' ...']; end
                                stradd{end+1}=sprintf('%s = %s', str1(n).name, regexprep(str2,'\s+',' '));
                            end
                        end
                    end
                end
            end
        end
    end
else
    if isempty(me), %Matlab<=2007a
        me=lasterror;
        str=regexprep(char(me.message),'[\t ]+',' ');
    else %2007b<=Matlab<=2008a
        str=regexprep(char(getReport(me)),'[\t ]+',' ');
    end
    if strfind(str,'<DisregardMessage>'), str='Process terminated by user'; PrimaryMessage=str;
    else
        idx0=strfind(str,'<nodetailsflag>');
        if ~isempty(idx0),
            str=str(idx0(1)+numel('<nodetailsflag>'):end);
        else
            str=regexprep(char(getReport(me)),'[\t ]+',' ');
        end
    end
end
errdatabase={...
    'error using ==> load','Potential missing data files or incomplete analyses. Try loading your project again (this will check for potential missing files) or re-running some of the analysis steps if some may have not completed correctly';
    'open image file','Potential missing data files or incomplete analyses. Try loading your project again (this will check for potential missing files) or re-running some of the analysis steps if some may have not completed correctly';
    'error using ==> fwrite|error using ==> fclose|error using ==> save','Try checking if storage space is full. Other potential issues include problems accesing a NFS drive, or incorrect file or folder permissions';
    'help memory','Insufficient memory, or memory fragmentation after long Matlab sessions. Try restarting matlab and attempting this procedure again';
    'index exceeds matrix dimensions','Potential missing Setup information. Try revising Setup fields for potential missing information';
    'file too small','Potential corrupt or incomplete volume. Try revising or regenerating source volumes';
    'opengl','Try using software OpenGL rendering (instead of hardware OpenGL). On unix machines you will need to issue the command ''opengl software'' right after starting Matlab (before any plotting functions are run)';
    };
str2=errdatabase(cellfun(@(a)~isempty(strfind(lower(str),a)),errdatabase(:,1)),2);
tb={};
if isdeployed, tb={'SPM compiled'};
else
    try, tb={spm('ver')}; tbt=spm('TBs'); if ~isempty(tbt), tb={[tb{:} ' +' sprintf(' %s',tbt.name)]}; end
    end
end
if iscell(connver)&&numel(connver)==2, str={'ERROR DESCRIPTION:',' ',str,['CONN',connver{2}],tb{:},['Matlab v.',version('-release')],['project: CONN',connver{1}]};
else str={'ERROR DESCRIPTION:',' ',str,['CONN v.',connver],tb{:},['Matlab v.',version('-release')]};
end

if ~isempty(projectname)
    try
        a=java.io.File(fileparts(projectname));
        k2=a.getTotalSpace;
        k0=a.canWrite;
        k=a.getUsableSpace;
        clear a;
        if k2==0, str=[str,{'storage: disconnected or unavailable'}];
        elseif k0, str=[str,{sprintf('storage: %.1fGb available',k*1e-9)}];
        else str=[str,{sprintf('storage: %.1fGb available (read-only)',k*1e-9)}];
        end
    end
end
try, str=[str {' '} regexp(evalc('conn_checkdistributionfiles([],''nolog'')'),'\n','split')]; end
if ~isempty(str2), str=[str,{' ','SUGGESTIONS:',' '},str2']; end
if ~isempty(stradd), str=[str {' ','stack details:',stradd{:}}]; end

end
