function [fh,TXT]=conn_reference(opts, varargin)
% CONN_REFERENCE automatic creation of "methods" section text with
% descriptions and citations of preprocessing/denoising/analysis 
% methods used
%
% conn_reference
% conn_reference('init');
% conn_reference('preprocessing');
% conn_reference('denoising')
% conn_reference('firstlevel');
% conn_reference('secondlevel');
%

global CONN_x
persistent connversions;

fields=struct('preproclog',[],...
              'denoisinginfo',[],...
              'firstlevelinfo',[],...
              'secondlevelinfo',[],...
              'fileout','referencesfile.html');

if nargin<1||isempty(opts), opts=''; end
for n=1:2:numel(varargin), assert(isfield(fields,varargin{n}),'unable to recognize option %s',varargin{n}); fields.(varargin{n})=varargin{n+1}; end


TXT={};CITATIONS={};
if isempty(opts), opts={'init','preprocessing','denoising','firstlevel'}; end
if ischar(opts), opts={opts}; end
for nopt=1:numel(opts),
    opt=opts{nopt};
    switch(lower(opt))
        case 'init'
            TTXT={};
            lines=conn_loadcfgfile(fullfile(fileparts(which(mfilename)),'conn_reference.txt'));
            connversions=lines;
            fnames=reshape(fieldnames(connversions),1,[]);
            fnames=sort(str2double(regexprep(fnames(cellfun('length',regexp(fnames,'VERSION\d+'))>0),'VERSION','')));
            for n=reshape(setdiff(1:max(fnames),fnames),1,[]), [nill,idx]=max((fnames<=n).*(1:numel(fnames))); connversions.(['VERSION',num2str(n)])=connversions.(['VERSION',num2str(fnames(idx))]); end
            if isfield(CONN_x,'ver'), conn_version=CONN_x.ver; 
            else conn_version=conn('ver');
            end
            lines.CITATION0=connversions.(['VERSION',regexprep(conn_version,'\..*','')]);
            options=struct;
            for line={'line_begin'};
                txt={}; citations={};
                if isfield(lines,line{1}), [txt,citations]=conn_reference_processtext(lines.(line{1}), lines, CITATIONS, options, conn_version);  end
                if ~isempty(txt), TTXT=[TTXT txt]; end
                if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
            end
            TTXT=conn_reference_singleparagraph(TTXT,false);
            TXT=[TXT, TTXT];

        case {'preproc','preprocessing'}
            TTXT={};
            if isempty(fields.preproclog),
                if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'log')||isempty(CONN_x.SetupPreproc.log), options=[];
                else
                    idx1=[];
                    if numel(CONN_x.SetupPreproc.log)>1,
                        str1={}; idx1=[];
                        for n=1:numel(CONN_x.SetupPreproc.log),
                            if isequal(CONN_x.SetupPreproc.log{n}{1},'timestamp'), str1{end+1}=CONN_x.SetupPreproc.log{n}{2}; idx1(end+1)=n; end
                        end
                        if numel(idx1)>1
                            answ=listdlg('Promptstring','Select preprocessing pipeline','selectionmode','single','liststring',str1,'ListSize',[320 300]);
                            if ~isempty(answ), idx1=idx1(answ); end
                        end
                    end
                    if isempty(idx1), idx1=numel(CONN_x.SetupPreproc.log); end
                    options=conn_reference_convertcell2struct(CONN_x.SetupPreproc.log{idx1});
                end
            else options=conn_reference_convertcell2struct(fields.preproclog);
            end
            if ~isempty(options)
                options.liststeps={}; % list steps in human-readable format
                [name1,name2]=conn_setup_preproc('steps');
                name2=regexprep(name2,'(\s*\(.*\).*$)|(^functional )|(^structural )','');
                name3=regexprep(lower(options.steps),'^run_|^update_|^interactive_','');
                for n=1:numel(name3)
                    if isempty(regexp(name3{n},'^functional_label|^functional_load|^functional_manual|^structural_manual|^functional_center|^structural_center|^structural_'))
                        [ok,idx]=ismember(name3{n},name1);
                        if ok, options.liststeps{end+1}=[lower(name2{idx}(1)) name2{idx}(2:end)]; end
                    end
                end
                if numel(options.liststeps)>1, options.liststeps=[' including ',sprintf('%s, ',options.liststeps{1:end-1}), ' and ',options.liststeps{end}];
                elseif numel(options.liststeps)==1, options.liststeps=[' including ',char(options.liststeps)];
                end

                lines=conn_loadcfgfile(fullfile(fileparts(which(mfilename)),'conn_reference_preproc.txt'));
                if isfield(options,'ver_CONN'), conn_version=options.ver_CONN;
                else conn_version=conn('ver');
                end
                lines.CITATION0=connversions.(['VERSION',regexprep(conn_version,'\..*','')]);
                for n1=double(numel(options.steps)<=1):numel(options.steps)
                    txt={}; citations={};
                    if n1==0
                        [txt,citations]=conn_reference_processtext(lines.default, lines, CITATIONS, options, conn_version);
                    else
                        stepname=regexprep(options.steps{n1},{'^run_','&'},{'','and'});
                        switch(stepname)
                            case fieldnames(lines)
                                [txt,citations]=conn_reference_processtext(lines.(stepname), lines, CITATIONS, options, conn_version);
                        end
                    end
                    if ~isempty(txt)&&isempty(TTXT), txt{1}=['Preprocessing: ',txt{1}]; end
                    if ~isempty(txt), TTXT=[TTXT txt]; end
                    if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                end
                TTXT=conn_reference_singleparagraph(TTXT);
                TXT=[TXT, TTXT];
            end

        case 'denoising'
            TTXT={};
            
            if isempty(fields.denoisinginfo)&&~CONN_x.isready(2), options=[];
            elseif isempty(fields.denoisinginfo), options=conn_reference_createdenoisingstruct(CONN_x.Preproc);
            else options=conn_reference_createdenoisingstruct(fields.denoisinginfo);
            end
            if ~isempty(options)
                lines=conn_loadcfgfile(fullfile(fileparts(which(mfilename)),'conn_reference_denoising.txt'));
                if isfield(CONN_x,'ver'), conn_version=CONN_x.ver;
                else conn_version=conn('ver');
                end
                lines.CITATION0=connversions.(['VERSION',regexprep(conn_version,'\..*','')]);
                if isinf(options.filter(2))&&options.filter(1)>0, bandpass={'highpass'}; options.filtered=1;
                elseif ~isinf(options.filter(2))&&options.filter(1)>0, bandpass={'bandpass'}; options.filtered=1;
                elseif ~isinf(options.filter(2))&&options.filter(1)==0, bandpass={'lowpass'}; options.filtered=1;
                else bandpass={}; options.filtered=0;
                end
                for line=[{'line_begin'},options.confounds_names,{'detrending',bandpass{:}}];
                    txt={}; citations={};
                    if isfield(lines,line{1}), [txt,citations]=conn_reference_processtext(lines.(line{1}), lines, CITATIONS, options, conn_version);
                    elseif isfield(lines,'otherwise')
                        tlines=lines;
                        tlines.(line{1})=regexprep(lines.otherwise,'otherwise',line{1});
                        [txt,citations]=conn_reference_processtext(tlines.(line{1}), tlines, CITATIONS, options, conn_version);
                    end
                    if ~isempty(txt)&&isempty(TTXT), txt{1}=['Denoising: ',txt{1}]; end
                    if ~isempty(txt), TTXT=[TTXT txt]; end
                    if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                end
                if isempty(bandpass)
                    TTXT=[TTXT(1), conn_reference_singleline(TTXT(2:end))];
                    TTXT={TTXT{1}, [TTXT{2},'.']};
                else
                    TTXT=[TTXT(1), conn_reference_singleline(TTXT(2:end-1)), TTXT(end)];
                    TTXT={TTXT{1}, [TTXT{2},', ',TTXT{3},'.']};
                end
                other={};
                if any(ismember({'wm','csf'},options.confounds_names)), other{end+1}='acompcor'; end
                if 1,
                    idx=find(cellfun('length',regexp(CONN_x.Setup.l2covariates.names,'^QC_DOF_session\d+$')));
                    if ~isempty(idx)
                        DOF2=zeros(1,CONN_x.Setup.nsubjects);
                        for nidx=1:numel(idx)
                            for nsub=1:CONN_x.Setup.nsubjects,
                                if ~isnan(CONN_x.Setup.l2covariates.values{nsub}{idx(nidx)}),
                                    DOF2(nsub)=DOF2(nsub)+CONN_x.Setup.l2covariates.values{nsub}{idx(nidx)};
                                end
                            end
                        end
                        options.QC_DOF_mean=round(10*mean(DOF2))/10;
                        options.QC_DOF_std=round(10*std(DOF2))/10;
                        options.QC_DOF_min=round(10*min(DOF2))/10;
                        options.QC_DOF_max=round(10*max(DOF2))/10;
                        other{end+1}='line_summary1';
                    end
                end
                for line=other
                    txt={}; citations={};
                    if isfield(lines,line{1}), [txt,citations]=conn_reference_processtext(lines.(line{1}), lines, CITATIONS, options, conn_version); end
                    if ~isempty(txt), TTXT=[TTXT txt]; end
                    if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                end
                TTXT=conn_reference_singleparagraph(TTXT,false);
                TXT=[TXT, TTXT];
            end

        case {'firstlevel','first-level'}
            TTXT={};
            lines=conn_loadcfgfile(fullfile(fileparts(which(mfilename)),'conn_reference_firstlevel.txt'));
            if isfield(CONN_x,'ver'), conn_version=CONN_x.ver; 
            else conn_version=conn('ver');
            end
            lines.CITATION0=connversions.(['VERSION',regexprep(conn_version,'\..*','')]);
            options=struct;
            options.multipleconditions=numel(CONN_x.Setup.conditions.names(1:end-1))>1;
            options.volumespace=CONN_x.Setup.spatialresolution<4;
            if isempty(fields.firstlevelinfo), 
                matchAnalyses=1:numel(CONN_x.Analyses);
                matchvvAnalyses=1:numel(CONN_x.vvAnalyses);
                matchdynAnalyses=1:numel(CONN_x.dynAnalyses);
            else
                matchAnalyses=find(ismember({CONN_x.Analyses.name},cellstr(fields.firstlevelinfo)));
                matchvvAnalyses=find(ismember({CONN_x.vvAnalyses.name},cellstr(fields.firstlevelinfo)));
                matchdynAnalyses=find(ismember({CONN_x.dynAnalyses.name},cellstr(fields.firstlevelinfo)));
            end
            for n=reshape(matchAnalyses,1,[])
                txt={}; citations={};
                toptions=CONN_x.Analyses(n);
                toptions.multipleconditions=options.multipleconditions;
                toptions.volumespace=options.volumespace;
                if ischar(CONN_x.Analyses(n).modulation),line={'TMOD'};      toptions.anyderivs=max([CONN_x.Analyses(n).regressors.deriv{:}]); if toptions.anyderivs>0, line{end+1}='ADDDERIV'; end; toptions.anyfbands=max([CONN_x.Analyses(n).regressors.fbands{:}]>1); if toptions.anyfbands>1, line{end+1}='ADDFBAND'; end; if options.multipleconditions, line{end+1}='ADDTMODWEIGHTS'; end
                elseif CONN_x.Analyses(n).modulation>0,  line={'gPPI'};      toptions.anyderivs=max([CONN_x.Analyses(n).regressors.deriv{:}]); if toptions.anyderivs>0, line{end+1}='ADDDERIV'; end; toptions.anyfbands=max([CONN_x.Analyses(n).regressors.fbands{:}]>1); if toptions.anyfbands>1, line{end+1}='ADDFBAND'; end; toptions.multipleconditions=1;
                elseif CONN_x.Analyses(n).type==1,       line={'RRC'};       toptions.anyderivs=max([CONN_x.Analyses(n).regressors.deriv{:}]); if toptions.anyderivs>0, line{end+1}='ADDDERIV'; end; toptions.anyfbands=max([CONN_x.Analyses(n).regressors.fbands{:}]>1); if toptions.anyfbands>1, line{end+1}='ADDFBAND'; end; line{end+1}='ADDWEIGHTS';
                elseif CONN_x.Analyses(n).type==2,       line={'SBC'};       toptions.anyderivs=max([CONN_x.Analyses(n).regressors.deriv{:}]); if toptions.anyderivs>0, line{end+1}='ADDDERIV'; end; toptions.anyfbands=max([CONN_x.Analyses(n).regressors.fbands{:}]>1); if toptions.anyfbands>1, line{end+1}='ADDFBAND'; end; line{end+1}='ADDWEIGHTS';
                elseif CONN_x.Analyses(n).type==3,       line={'SBCRRC'};    toptions.anyderivs=max([CONN_x.Analyses(n).regressors.deriv{:}]); if toptions.anyderivs>0, line{end+1}='ADDDERIV'; end; toptions.anyfbands=max([CONN_x.Analyses(n).regressors.fbands{:}]>1); if toptions.anyfbands>1, line{end+1}='ADDFBAND'; end; line{end+1}='ADDWEIGHTS';
                end
                if isfield(toptions,'sources')&&~isempty(toptions.sources) % finished analyses only
                    if ~options.multipleconditions, toptions.weight=-toptions.weight; end % marks single-condition analyses
                    if numel(toptions.sources)==1, toptions.listseeds=toptions.sources{1};
                    elseif numel(toptions.sources)<=4, toptions.listseeds=[toptions.sources{1},sprintf(', %s',toptions.sources{2:end-1}),' and ',toptions.sources{end}];
                    elseif all(cellfun('length',regexp(toptions.sources,'^networks\.'))), toptions.listseeds=sprintf('%d HPC-ICA network ROIs [CITATION0]',numel(toptions.sources));
                    elseif all(cellfun('length',regexp(toptions.sources,'^atlas\.'))), toptions.listseeds=sprintf('%d Harvard-Oxford atlas ROIs [CITATION2]',numel(toptions.sources));
                    elseif all(cellfun('length',regexp(toptions.sources,'^atlas\.|^networks\.'))), toptions.listseeds=sprintf('%d HPC-ICA networks [CITATION0] and Harvard-Oxford atlas ROIs [CITATION2]',numel(toptions.sources));
                    else toptions.listseeds=sprintf('%d ROIs',numel(toptions.sources));
                    end
                    if numel(toptions.conditions)==0, toptions.listconditions='no task or experimental';
                    elseif numel(toptions.conditions)==1, toptions.listconditions=toptions.conditions{1};
                    elseif numel(toptions.conditions)<=4, toptions.listconditions=[toptions.conditions{1},sprintf(', %s',toptions.conditions{2:end-1}),' and ',toptions.conditions{end}];
                    else toptions.listconditions=sprintf('%d task or experimental',numel(toptions.conditions));
                    end
                    LTXT={};
                    for nline=1:numel(line)
                        if isfield(lines,line{nline})
                            [txt,citations]=conn_reference_processtext(lines.(line{nline}), lines, CITATIONS, toptions, conn_version);
                            if ~isempty(txt), LTXT=[LTXT txt]; end
                            if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                        end
                    end
                    if ~isempty(LTXT), TTXT=[TTXT conn_reference_singleparagraph(LTXT,false)]; end
                end
            end
            
            for n=reshape(matchvvAnalyses,1,[])
                switch(CONN_x.vvAnalyses(n).regressors.measuretype{1})
                    case 1, 
                        if CONN_x.vvAnalyses(n).regressors.global{1}==0&&CONN_x.vvAnalyses(n).regressors.deriv{1}==0, line={'LCOR'};
                        elseif CONN_x.vvAnalyses(n).regressors.global{1}>0&&CONN_x.vvAnalyses(n).regressors.deriv{1}==0, line={'IC'};
                        elseif CONN_x.vvAnalyses(n).regressors.global{1}==0&&CONN_x.vvAnalyses(n).regressors.deriv{1}>0, line={'RCOR'};
                        elseif CONN_x.vvAnalyses(n).regressors.global{1}>0&&CONN_x.vvAnalyses(n).regressors.deriv{1}>0, line={'RSIM'};
                        end
                    case 5, line={'GCOR'};
                    case 2, line={'fcMVPA'};
                    case 3, line={'groupICA'};
                    case 4, line={'groupPCA'};
                    case 6, line={'ALFF'};
                    case 7, line={'fALFF'};
                    case 8, line={'IHC'};
                end
                txt={}; citations={};
                toptions=CONN_x.vvAnalyses(n);
                toptions.multipleconditions=options.multipleconditions;
                toptions.volumespace=options.volumespace;
                for tfield=reshape(fieldnames(CONN_x.vvAnalyses(n).regressors),1,[]),
                    toptions.(tfield{1})=CONN_x.vvAnalyses(n).regressors.(tfield{1}){1};
                end
                toptions.finitedimensions=~isinf(toptions.dimensions_in);
                toptions.ismasked=~isempty(toptions.mask);
                toptions.bpfilter=CONN_x.Preproc.filter;
                if toptions.measuretype==3||toptions.measuretype==4
                    if isempty(toptions.options), toptions.options='GICA3tanh'; end
                    toptions.algtype=regexprep(regexprep(toptions.options,'GICA\d',''),{'^tanh$','^pow3$','^gauss$'},{'hyperbolic tangent (G1)','cubic (G3)','gaussian (G2)'});
                    toptions.bprtype=regexp(toptions.options,'GICA\d','match','once');
                end
                if 1,
                    for nline=1:numel(line)
                        if isfield(lines,line{nline})
                            [txt,citations]=conn_reference_processtext(lines.(line{nline}), lines, CITATIONS, toptions, conn_version);
                        end
                    end
                    if ~isempty(txt), TTXT=[TTXT txt]; end
                    if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                end
            end

            for n=reshape(matchdynAnalyses,1,[])
                line={'dynICA'};
                txt={}; citations={};
                toptions=CONN_x.dynAnalyses(n);
                toptions.multipleconditions=options.multipleconditions;
                toptions.volumespace=options.volumespace;
                for tfield=reshape(fieldnames(CONN_x.dynAnalyses(n).regressors),1,[]),
                    toptions.(tfield{1})=CONN_x.dynAnalyses(n).regressors.(tfield{1}){1};
                end
                if numel(toptions.sources)==1, toptions.listseeds=toptions.sources{1};
                elseif numel(toptions.sources)<=4, toptions.listseeds=[toptions.sources{1},sprintf(', %s',toptions.sources{2:end-1}),' and ',toptions.sources{end}];
                elseif all(cellfun('length',regexp(toptions.sources,'^networks\.'))), toptions.listseeds=sprintf('%d HPC-ICA network ROIs [CITATION0]',numel(toptions.sources));
                elseif all(cellfun('length',regexp(toptions.sources,'^atlas\.'))), toptions.listseeds=sprintf('%d Harvard-Oxford atlas ROIs [CITATION2]',numel(toptions.sources));
                elseif all(cellfun('length',regexp(toptions.sources,'^atlas\.|^networks\.'))), toptions.listseeds=sprintf('%d HPC-ICA networks [CITATION0] and Harvard-Oxford atlas ROIs [CITATION2]',numel(toptions.sources));
                else toptions.listseeds=sprintf('%d ROIs',numel(toptions.sources));
                end
                if 1,
                    for nline=1:numel(line)
                        if isfield(lines,line{nline})
                            [txt,citations]=conn_reference_processtext(lines.(line{nline}), lines, CITATIONS, toptions, conn_version);
                        end
                    end
                    if ~isempty(txt), TTXT=[TTXT txt]; end
                    if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                end
            end

            %TTXT=conn_reference_singleparagraph(TTXT,false);
            TXT=[TXT, TTXT];

        case {'secondlevel','second-level'}
            TTXT={};
            lines=conn_loadcfgfile(fullfile(fileparts(which(mfilename)),'conn_reference_secondlevel.txt'));
            if isfield(CONN_x,'ver'), conn_version=CONN_x.ver; 
            else conn_version=conn('ver');
            end
            lines.CITATION0=connversions.(['VERSION',regexprep(conn_version,'\..*','')]);
            options=struct;
            options.multipleconditions=numel(CONN_x.Setup.conditions.names(1:end-1))>1;
            options.volumespace=CONN_x.Setup.spatialresolution<4;
            for n=0:numel(fields.secondlevelinfo)
                txt={}; citations={};
                toptions=options;
                if n==0, line={'init'};
                else line=fields.secondlevelinfo(n);
                end
                if (n==0&&any(ismember(fields.secondlevelinfo,{'RFT','RANDOM','TFCE'})))||(n>0&&ismember(fields.secondlevelinfo{n},{'RFT','RANDOM','TFCE'})), toptions.voxel='voxel';toptions.Voxel='Voxel';
                else toptions.voxel='connection';toptions.Voxel='Connection';
                end
                LTXT={};
                for nline=1:numel(line)
                    if isfield(lines,line{nline})
                        [txt,citations]=conn_reference_processtext(lines.(line{nline}), lines, CITATIONS, toptions, conn_version);
                        if ~isempty(txt), LTXT=[LTXT txt]; end
                        if ~isempty(citations), CITATIONS=[CITATIONS citations]; end
                    end
                end
                if ~isempty(LTXT), TTXT=[TTXT conn_reference_singleparagraph(LTXT,false)]; end
            end
            TTXT=conn_reference_singleparagraph(TTXT,false);
            TXT=[TXT, TTXT];

    end
end
if ~isempty(TXT)
    txt=[TXT, {''}, CITATIONS];
    if ~nargout, fprintf('%s\n',txt{:}); end
    if isempty(fields.fileout), 
        fh={};
        fh{end+1}=sprintf('<html><body><p style="font-size:16px"><b>Methods</b></p>\n');
        for n=1:numel(TXT), 
            fh{end+1}=sprintf('<p style="font-size:12px;font-family:helvetica;line-height:2">\n');
            txt=TXT{n};
            txt=regexprep(txt,{'\&','<','>'},{'&amp;','&lt;','&gt;'});
            txt=regexprep(txt,'\s*\[(\d+)\]','<sup>[$1]</sup>');
            txt=regexprep(txt,'\]<\/sup><sup>\[',',');
            txt=regexprep(txt,'^Preprocessing: ','<b>Preprocessing</b>: ');
            txt=regexprep(txt,'^Denoising: ','<b>Denoising</b>: ');
            txt=regexprep(txt,'^First-level analysis (.*)?\: ','<b>First-level analysis</b> $1: ');
            txt=regexprep(txt,'^Group-level analyses','<b>Group-level analyses</b>');
            txt=regexprep(txt,'\n','');
            %txt=regexprep(txt,'CONN|SPM|ART','<em>$0</em>');
            fh{end+1}=sprintf('%s</p>\n',txt);
        end
        fh{end+1}=sprintf('<p style="font-size:16px"><b>References</b></p>\n');
        for n=1:numel(CITATIONS), 
            txt=CITATIONS{n};
        fh{end+1}=sprintf('');
            txt=regexprep(txt,{'\&','<','>'},{'&amp;','&lt;','&gt;'});
            txt=regexprep(txt,'\s*\[(\d+)\]','<sup>[$1]</sup>');
            fh{end+1}=sprintf('<p style="font-size:11px;font-family:helvetica;line-height:2">%s</p>\n',txt);
        end
        fh{end+1}=sprintf('</body></html>\n');
        fh=[fh{:}];
    else
        fh=fopen(fields.fileout,'wt');
        fprintf(fh,'<html><body><p style="font-size:12px;font-family:helvetica;line-height:2"><center>Copy and paste the section below to your manuscript <b>Methods</b> section. This text is distributed under a Public Domain Dedication license (<a href="https://creativecommons.org/publicdomain/zero/1.0/">CC0 1.0</a>) and <a href="https://web.conn-toolbox.org/resources/citing-conn">it can be used, copied, modified, and distributed freely</a>.</center></p><br></br>\n');
        fprintf(fh,'<p style="font-size:16px"><b>Methods</b></p>\n');
        for n=1:numel(TXT), 
            fprintf(fh,'<p style="font-size:12px;font-family:helvetica;line-height:2">\n');
            txt=TXT{n};
            txt=regexprep(txt,{'\&','<','>'},{'&amp;','&lt;','&gt;'});
            txt=regexprep(txt,'\s*\[(\d+)\]','<sup>[$1]</sup>');
            txt=regexprep(txt,'\]<\/sup><sup>\[',',');
            txt=regexprep(txt,'^Preprocessing: ','<b>Preprocessing</b>: ');
            txt=regexprep(txt,'^Denoising: ','<b>Denoising</b>: ');
            txt=regexprep(txt,'^First-level analysis (.*)?\: ','<b>First-level analysis</b> $1: ');
            txt=regexprep(txt,'^Group-level analyses','<b>Group-level analyses</b>');
            txt=regexprep(txt,'\n','');
            %txt=regexprep(txt,'CONN|SPM|ART','<em>$0</em>');
            fprintf(fh,'%s</p>\n',txt);
        end
        fprintf(fh,'<p style="font-size:16px"><b>References</b></p>\n');
        for n=1:numel(CITATIONS), 
            txt=CITATIONS{n};
            txt=regexprep(txt,{'\&','<','>'},{'&amp;','&lt;','&gt;'});
            txt=regexprep(txt,'\s*\[(\d+)\]','<sup>[$1]</sup>');
            fprintf(fh,'<p style="font-size:11px;font-family:helvetica;line-height:2">%s</p>\n',txt);
        end
        fprintf(fh,'</body></html>\n');
        fclose(fh);
        fh=conn_fullfile(fields.fileout);
    end
end
end

function b=conn_reference_convertcell2struct(a)
for n=1:2:numel(a)-1
    b.(a{n})=a{n+1};
end
end

function b=conn_reference_createdenoisingstruct(a)
b=a;
b.confounds_names={};
taskeffects=sum(cellfun('length',regexp(a.confounds.names,'^Effect of '))>0);
for n=1:numel(a.confounds.names)
    fname='';
    if n==0
    else
        switch(a.confounds.names{n})
            case 'Grey Matter', fname='gm';
            case 'White Matter', fname='wm';
            case 'CSF', fname='csf';
            otherwise
                if taskeffects==1&&~isempty(regexp(a.confounds.names{n},'^Effect of rest$')), fname='session';
                elseif ~isempty(regexp(a.confounds.names{n},'^Effect of ')), fname='task';
                else fname=regexprep(a.confounds.names{n},'[^\w\d]','_');
                end
        end
        if ~ismember(fname,b.confounds_names)
            b.confounds_names{end+1}=fname;
            b.(['confounds_',fname,'_','deriv'])=a.confounds.deriv{n};
            b.(['confounds_',fname,'_','power'])=a.confounds.power{n};
            b.(['confounds_',fname,'_','dimensions'])=min(a.confounds.dimensions{n})*a.confounds.power{n}*(1+a.confounds.deriv{n});
            b.(['confounds_',fname,'_','filter'])=a.confounds.filter{n};
        else
            b.(['confounds_',fname,'_','deriv'])=max(b.(['confounds_',fname,'_','deriv']),a.confounds.deriv{n});
            b.(['confounds_',fname,'_','power'])=max(b.(['confounds_',fname,'_','power']),a.confounds.power{n});
            b.(['confounds_',fname,'_','dimensions'])=b.(['confounds_',fname,'_','dimensions'])+min(a.confounds.dimensions{n})*a.confounds.power{n}*(1+a.confounds.deriv{n});
            b.(['confounds_',fname,'_','filter'])=max(b.(['confounds_',fname,'_','filter']),a.confounds.filter{n});
        end
    end
end
end

function txtout=conn_reference_singleparagraph(txtin,dolast)
if nargin<2, dolast=true; end
txtout='';
for n=1:numel(txtin),
    if n==1, txtout=[txtin{n}]; 
    elseif n<numel(txtin), txtout=[txtout, ' ', txtin{n}]; 
    elseif dolast&&numel(txtin)>2, txtout=[txtout, ' Last, ',lower(txtin{n}(1)),txtin{n}(2:end)];
    else txtout=[txtout, ' ',txtin{n}];
    end
end
txtout={txtout};
end

function txtout=conn_reference_singleline(txtin)
txtout='';
for n=1:numel(txtin),
    if n==1, txtout=[txtin{n}]; 
    elseif n<numel(txtin), txtout=[txtout, ', ', txtin{n}]; 
    else txtout=[txtout, ', and ',txtin{n}];
    end
end
txtout={txtout};
end

function [line,citation]=conn_reference_processtext(line, lines, citations, options, conn_version)
if nargin<5||isempty(conn_version), conn_version=conn('ver'); end
if ~iscell(line), line={line}; end
citation={};
for nline=1:numel(line)
    ok=false;
    while ~ok
        key=regexp(line{nline},'\[[^\[\]]*\]','match','once');
        if isempty(key), ok=true;
        else
            keyout='';
            switch(key)
                case '[SPMVERSION]', keyout=spm('ver');
                case '[SPMFULLVERSION]', keyout=spm('version');
                case '[SPMRELEASE]', keyout=regexprep(spm('version'),'\D+(\d+)\s*\((.*?)\)','$1.$2');
                case '[CONNVERSION]', keyout=conn_version;
                case '[CONNFULLVERSION]', keyout=sprintf('CONN (%s)',conn_version);
                case '[CONNRELEASE]', keyout=conn_version;
                case '[VALUE 1 ver_SPM]', if isfield(options,'ver_SPM'), keyout=regexprep(options.ver_SPM,'\s.*$',''); else keyout=spm('ver'); end
                otherwise
                    if ~isempty(regexp(key,'^\[CITATION\d+')) % adds a citation
                        txt=char(lines.(regexprep(key,'\[|\]','')));
                        inlist=ismember(regexprep(citations,'^[{\[]\d+[}\]] ',''),txt);
                        if any(inlist)
                            inlist=find(inlist,1);
                            keyout=sprintf('{%d}',str2double(regexprep(citations{inlist},'^[{\[](\d+)[}\]].*$','$1')));
                        else
                            keyout=sprintf('{%d}',numel(citations)+1);
                            citation{end+1}=[sprintf('{%d} ',numel(citations)+1) txt];
                            citations{end+1}=citation{end};
                        end
                    elseif ~isempty(regexp(key,'^\[VALUE\s\d+ ')) % adds a value specified in the log
                        fidx=regexprep(key,'^\[VALUE\s(\d+)\s.*\]$','$1');
                        fname=regexprep(key,'^\[VALUE\s\d+\s(.*)\]$','$1');
                        fvalue=options.(fname);
                        if ~isempty(fvalue)&&~ischar(fvalue), fvalue=fvalue(str2double(fidx)); end
                        if ischar(fvalue), keyout=fvalue;
                        elseif iscell(fvalue)&&numel(fvalue)==1&&isnumeric(fvalue{1}), keyout=mat2str(fvalue{1});
                        elseif iscell(fvalue), keyout=char(fvalue);
                        else keyout=mat2str(fvalue);
                        end
                    elseif ~isempty(regexp(key,'^\[IF ')) % adds conditional to a value in the log
                        fname=regexprep(key,'^\[IF\s([^\s]+)\s.*\]$','$1');
                        fvalue=regexprep(key,'^\[IF\s[^\s]+\s([^\s]+)\s.*\]$','$1');
                        ftext= regexprep(key,'^\[IF\s[^\s]+\s[^\s]+\s(.*)\]$','$1');
                        if ~isnan(str2double(fvalue)), dokeyout=any(options.(fname)==str2double(fvalue));
                        else, dokeyout=any(ismember(options.(fname),fvalue));
                        end
                        if dokeyout, keyout=ftext; end
                    elseif ~isempty(regexp(key,'^\[IFEMPTY ')) % adds conditional to a value in the log
                        fname=regexprep(key,'^\[IFEMPTY\s([^\s]+)\s.*\]$','$1');
                        fvalue=regexprep(key,'^\[IFEMPTY\s[^\s]+\s([^\s]+)\s.*\]$','$1');
                        ftext= regexprep(key,'^\[IFEMPTY\s[^\s]+\s[^\s]+\s(.*)\]$','$1');
                        dokeyout=isempty(options.(fname))==str2double(fvalue);
                        if dokeyout, keyout=ftext; end
                    else fprintf('warning: unrecognized key %s\n',key);
                    end
            end
            line{nline}=regexprep(line{nline},'\[[^\[\]]*\]',keyout,'once');
        end
    end
end
line=regexprep(line,'\{(\d+)\}','[$1]');
citation=regexprep(citation,'\{(\d+)\}','[$1]');
end
