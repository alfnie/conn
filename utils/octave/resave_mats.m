CONNDIR = fileparts(which('conn'));

for fn = cellstr(spm_select('FPListRec',CONNDIR,'.*\.mat$'))'
    fprintf('Processing %s...',fn{1});
    vars = load(fn{1});
    nameVars = fieldnames(vars);
    cellfun(@(v) assignin('base',v,vars.(v)), nameVars);
    save(fn{1},nameVars{:},'-v7')
    fprintf('Done!\n')
end
