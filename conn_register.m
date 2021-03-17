function askagain=conn_register(option)

if nargin<1, option='register'; end
file_registered=fullfile(fileparts(which(mfilename)),'conn_register.txt');
if strcmp(option,'register')&&conn_existfile(file_registered)
    answ=conn_questdlg({'We would like to collect some simple demographic information from our users.','This serves to keep track of the toolbox impact and it helps in our funding applications.','Would you like to register now? (registration takes less than a minute)'},'CONN toolbox registration','Yes','No, ask me later','No, never ask again','No, ask me later');
    if strcmp(answ,'Yes'), option='forceregister'; end
    if strcmp(answ,'No, never ask again'), option='donotaskagain'; end
end
if strcmp(option,'forceregister'),
    web('http://www.alfnie.com/software/conn-toolbox-registration','-browser');
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
