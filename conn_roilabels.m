function [ROInames,ROIidx]=conn_roilabels(filename) 
% internal: reads ROI labels from atlas file

[roi_path_dir,roi_path_name,roi_path_ext]=fileparts(filename);
ok=false;
if conn_existfile(fullfile(roi_path_dir,[roi_path_name,'.txt'])),
    lines=conn_fileutils('textread',fullfile(roi_path_dir,[roi_path_name,'.txt']),'%s','delimiter','\n');
    for opts=1:3
        try
            switch(opts)
                case 1, % (ROI_NUMBER ROI_LABEL)
                    tlines=lines(cellfun('length',lines)>0);
                    words=regexp(tlines,'^\s*(\d+)\s+(.+)\s*$','tokens','once');
                    words=cat(1,words{:});
                    ROIidx=str2double(words(:,1));
                    ROInames=words(:,2);
                    assert(size(words,1)==numel(tlines));
                    assert(all(~isnan(ROIidx)));
                    assert(all(cellfun('length',ROInames)>0));
                    ok=true;
                case 2, % FreeSurfer *LUT.txt or equivalent format (ROI_NUMBER ROI_LABEL)
                    tlines=lines(cellfun('length',lines)>0);
                    tlines=regexprep(tlines,'^#.*','');
                    tlines=tlines(cellfun('length',tlines)>0);
                    words=regexp(tlines,'\s\t','split');
                    assert(all(cellfun('length',words)>=2));
                    words=cellfun(@(x)x(1:2),words,'uni',0);
                    words=cat(1,words{:});
                    ROIidx=str2double(words(:,1));
                    ROInames=words(:,2);
                    assert(all(~isnan(ROIidx)));
                    assert(all(cellfun('length',ROInames)>0));
                    ok=true;
                case 3, % (ROI_LABEL)
                    ROInames=lines;
                    ROIidx=1:numel(ROInames);
                    valid=cellfun('length',ROInames)>0;
                    ROInames=ROInames(valid);
                    ROIidx=ROIidx(valid);
                    ok=true;
            end
        end
        if ok, break; end
    end
elseif conn_existfile(fullfile(roi_path_dir,[roi_path_name,'.csv'])),
    lines=conn_fileutils('textread',fullfile(roi_path_dir,[roi_path_name,'.csv']),'%s','delimiter','\n');
    lines=lines(cellfun('length',lines)>0);
    for opts=1:4
        try
            switch(opts)
                case 1, % (ROI_LABEL,ROI_NUMBER) with headerline
                    tlines=lines(2:end);
                    words=regexp(tlines,'^\s*(.+)\s*,\s*(\d+)\s*$','tokens','once');
                    words=cat(1,words{:});
                    ROIidx=str2double(words(:,2));
                    ROInames=words(:,1);
                    assert(all(~isnan(ROIidx)));
                    assert(size(words,1)==numel(tlines));
                    ok=true;
                case 2, % (ROI_NUMBER,ROI_LABEL) with headerline
                    tlines=lines(2:end);
                    words=regexp(tlines,'^\s*(\d+)\s*,\s*(.+)\s*$','tokens','once');
                    words=cat(1,words{:});
                    ROIidx=str2double(words(:,1));
                    ROInames=words(:,2);
                    assert(all(~isnan(ROIidx)));
                    assert(size(words,1)==numel(tlines));
                    ok=true;
                case 3, % (ROI_LABEL,ROI_NUMBER) without headerline
                    tlines=lines;
                    words=regexp(tlines,'^\s*(.+)\s*,\s*(\d+)\s*$','tokens','once');
                    words=cat(1,words{:});
                    ROIidx=str2double(words(:,2));
                    ROInames=words(:,1);
                    assert(all(~isnan(ROIidx)));
                    assert(size(words,1)==numel(tlines));
                    ok=true;
                case 4, % (ROI_NUMBER,ROI_LABEL) without headerline
                    tlines=lines;
                    words=regexp(tlines,'^\s*(\d+)\s*,\s*(.+)\s*$','tokens','once');
                    words=cat(1,words{:});
                    ROIidx=str2double(words(:,1));
                    ROInames=words(:,2);
                    assert(all(~isnan(ROIidx)));
                    assert(size(words,1)==numel(tlines));
                    ok=true;
            end
        end
        if ok, break; end
    end
elseif conn_existfile(fullfile(roi_path_dir,[roi_path_name,'.xls'])),
    [idxpairs,PU]=xlsread(fullfile(roi_path_dir,[roi_path_name,'.xls'])); % ROI_NUMBER,ROI_LABEL
    if size(PU,1)==size(idxpairs,1), ROInames=PU(:,1);ROIidx=idxpairs(:,1); ok=true;
    elseif size(PU,1)==size(idxpairs,1)+1, ROInames=PU(2:end,1);ROIidx=idxpairs(:,1); ok=true;
    end
end
    
if ~ok
    ROInames={};
    ROIidx=[];
end
