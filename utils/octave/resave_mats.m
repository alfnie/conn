CONNDIR = fileparts(which('conn'));

for fn = dir(fullfile(CONNDIR,'**/*.mat'))'
    fnfull = fullfile(fn.folder, fn.name);
    fprintf('Processing %s...',fnfull);
    
    % Load content
    vars = load(fnfull);
    nameVars = fieldnames(vars);
    cellfun(@(v) assignin('base',v,vars.(v)), nameVars);
    
    % Re-save
    save(fnfull,nameVars{:},'-v7')
    
    % Clear
    cellfun(@clear, nameVars);

    fprintf('Done!\n')
end
