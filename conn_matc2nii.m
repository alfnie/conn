function [outfiles,outvals]=conn_matc2nii(filename,dowaitbar,localfolder,dofast)
% CONN_MATC2NII
% converts .mat/.matc volumes to .nii format

global CONN_x;
if nargin<2||isempty(dowaitbar), dowaitbar=1;end
if nargin<3||isempty(localfolder),localfolder='';end 
if nargin<4||isempty(dofast), dofast=1;end

if nargin<1||(isempty(filename)&&~iscell(filename)),
    filename={};
	filepathresults=CONN_x.folders.preprocessing;
	nconditions=length(CONN_x.Setup.conditions.names)-1;
    for nsub=1:CONN_x.Setup.nsubjects,
        for ncondition=1:nconditions,
            [icondition,isnew]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition});
            if isnew, error(['Mismatched condition ',CONN_x.Setup.conditions.names{ncondition}]); end
            filename{end+1}=fullfile(filepathresults,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition,'%03d'),'.mat']);
        end
    end
end
if ~iscell(filename), filename={filename};end
if any(conn_server('util_isremotefile',filename)), 
    [outfiles,outvals]=conn_server('run',mfilename,conn_server('util_localfile',filename),dowaitbar,localfolder,dofast); 
    outfiles=conn_server('util_remotefile',outfiles);
    return
end
filename=conn_server('util_localfile',filename);
filename=regexprep(filename,'\.matc$','.mat');

outfiles={}; outvals={};
if dowaitbar, warning off; h=conn_waitbar(0,['Converting to nifti. Please wait...']); warning on; end
for nfile=1:length(filename),
    Y=conn_vol(filename{nfile});
    if isfield(Y,'softlink')&&~isempty(Y.softlink),
        str1=regexp(Y.fname,'Subject\d+','match'); if ~isempty(str1), Y.softlink=regexprep(Y.softlink,'Subject\d+',str1{end}); end
        [file_path,file_name,file_ext]=fileparts(Y.fname);
        filename{nfile}=fullfile(file_path,Y.softlink);
    end
    if ~ismember(filename{nfile},outfiles)
        [filenamepath,filenamename,filenameext]=fileparts(filename{nfile});
        if isempty(localfolder), fname=conn_prepend('nifti',fullfile(filenamepath,[filenamename,'.nii']));
        else fname=conn_prepend('nifti',fullfile(localfolder,[filenamename,'.nii'])); spm_unlink(fname);
        end
        V=struct('fname',fname,...
            'mat',Y.matdim.mat,...
            'dim',[Y.matdim.dim],...
            'n',[1,1],...
            'pinfo',[1;0;0],...
            'dt',[spm_type('float32'),spm_platform('bigend')],...
            'descrip',mfilename);
        V=repmat(V,[Y.size.Nt,1]);for nt=1:Y.size.Nt,V(nt).n=[nt,1];end
        V=spm_create_vol(V);
        mV=zeros(1,Y.size.Nt);nV=0;
        if dofast
            for nslice=1:Y.size.Ns
                [X,idx]=conn_get_slice(Y,nslice);
                mV=mV+sum(X,2)';
                nV=nV+numel(idx);
                for nt=1:Y.size.Nt
                    Z=zeros(Y.matdim.dim(1:2));
                    Z(idx)=X(nt,:);
                    V(nt) = spm_write_plane(V(nt),Z,nslice);
                end
                if dowaitbar, conn_waitbar((nfile-1+nslice/Y.size.Ns)/length(filename),h); end
            end
            mV=mV/max(eps,nV);
        else
            for nt=1:Y.size.Nt
                [Z,idx]=conn_get_time(Y,nt);
                V(nt)=spm_write_vol(V(nt),Z);
                mV(nt)=mean(Z(idx));
            end
        end
        if ~isempty(localfolder), 
            fname2=conn_prepend('nifti',fullfile(filenamepath,[filenamename,'.nii']));
            fprintf('moving file from %s to %s...',fname,fname2);
            if ispc, [ok,nill]=system(['move "',fname,'" "',fname2,'"']);
            else, [ok,nill]=system(['mv -f ''',fname,''' ''',fname2,'''']);
            end
            if ~isequal(ok,0), error('Error moving file %s to %s, check target permissions',fname,fname2); end
            fprintf('done\n');
        end
        outfiles{end+1}=filename{nfile};
        outvals{end+1}=mV;
    end
    if dowaitbar&&~dofast, conn_waitbar(nfile/length(filename),h); end
end
if dowaitbar, close(h); end

