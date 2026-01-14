function askagain=conn_register(option)

if nargin<1, option='register'; end
file_registered=fullfile(fileparts(which(mfilename)),'conn_register.txt');
if strcmp(option,'register')&&conn_existfile(file_registered)
    answ=conn_questdlg({'Register to help us quantify CONN''s usage in the research community','and to receive CONN''s newsletter, with monthly news and updates about CONN.', ' ', 'Would you like to register now?',' ','[You may use CONN without registering. Registration is optional, free, and takes under a minute.]'},'CONN toolbox registration','Yes, register now','Not now','No, never ask again','Not now');
    if strcmp(answ,'Yes, register now'), option='forceregister'; end
    if strcmp(answ,'No, never ask again'), option='donotaskagain'; end
end
if strcmp(option,'forceregister'),
    web('https://www.conn-toolbox.org/registration','-browser');
    option='donotaskagain';
end
if strcmp(option,'donotaskagain')
    spm_unlink(file_registered);
    if conn_existfile(file_registered)
        conn_disp(['Warning: Unable to delete file ',file_registered]);
        conn_disp('Please delete this file manually to avoid being asked to register again');
    end
end
askagain=conn_existfile(file_registered);
