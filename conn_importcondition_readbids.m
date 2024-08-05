function [onsets,durations,conditions]=conn_importcondition_readbids(filename)
% [onsets,durations,conditions]=conn_importcondition_readbids(filename)
%   reads BIDS *_events.tsv file:
%     .tsv file (tab-separated fields)
%     one header line, and three columns: onset, duration, trial_type
%     (notes: any additional columns are disregarded; columns may appear in any order; "trial_type" column alternative names: "event_type", "condition")
%

str=fileread(filename);
str=regexp(str,'[\n\r]','split');
str=str(cellfun('length',str)>0);
str=regexp(str,'\t','split');
header=regexprep(str{1},'^\s+|\s+$',''); str=str(2:end);
nstr=cellfun('length',str);
idx1=find(strcmp(lower(header),'onset'),1);
if isempty(idx1), error('no ''onset'' field found in %s',filename); end
idx2=find(strcmp(lower(header),'duration'),1);
if isempty(idx2), error('no ''duration'' field found in %s',filename); end
idx3=find(strcmp(lower(header),'trial_type'),1);
if isempty(idx3), idx3=find(strcmp(lower(header),'event_type'),1); end
if isempty(idx3), idx3=find(strcmp(lower(header),'condition'),1); end
if isempty(idx3), error('no ''trial_type'' field found in %s',filename); end
%str=cat(1,str{nstr==numel(header)});
str=cat(1,str{nstr>=max([idx1 idx2 idx3])});

onsets=str2double(str(:,idx1));
durations=str2double(str(:,idx2));
conditions=str(:,idx3);
end

