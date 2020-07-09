function changed=conn_setpermissions(filemask_options)
% CONN_SETPERMISSIONS checks/sets file mode creation mask (umask)
%
% conn_setpermissions
%   checks/sets file mode creation mask (umask) to grant group-level read&write access (i.e. mask='002')
%
% conn_setpermissions(mask) 
%   sets the file mode creation mask to a value equal or smaller (bitwise) than mask
%   e.g. conn_setpermissions('002')
%

changed=false;
if nargin<1||isempty(filemask_options), filemask_options='0002'; end
try
    if isunix&&~ismac
        [ok,msg]=system('umask');
        if ~ok,
            msg1=base2dec(regexprep(msg,'[^012345678]',''),8);
            msg2=base2dec(filemask_options,8);
            if bitand(msg1,intmax-msg2)
                str1=dec2base(msg1,8);
                str2=dec2base(bitand(msg1,msg2),8);
                fprintf(...
                    'warning: unexpected file mode creation mask found (%s). Mask set to %s. Change in your user configuration files the string "umask %s" to "umask %s" in order to stop seeing this warning in the future',...
                    str1,str2,str1,str2);
                try
                    [ok,tmsg]=system('grep "umask " ~/.*');
                    tmsg=regexp(tmsg,'\n','split');
                    tmsg=tmsg(cellfun('length',tmsg)>0);
                    tmsg=tmsg(cellfun('length',regexp(tmsg,['umask ',regexprep(str1,'^0+','0*')]))>0);
                    if ~isempty(tmsg)
                        tmsg2=regexprep(tmsg,regexprep(str1,'^0+','0*'),str2);
                        fprintf(' e.g. change "%s" to "%s"',sprintf('%s ',tmsg{:}),sprintf('%s',tmsg2{:}));
                    end
                end
                fprintf('\n');
                [nill,msg]=system(sprintf('umask %s',str2));
                changed=true;
            end
        end
    end
end
