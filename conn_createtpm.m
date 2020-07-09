function newTPMfile=conn_createtpm(lesionmask, structural, newTPMfile)
% conn_createtpm creates new TPM.nii file by adding one new tissue class defined from a lesion mask
% 
% conn_createtpm(lesionmask, structural, newTPMfile)
%   lesionmask: input lesion mask file
%   structural: input reference structural file in same space as lesionmask
%   newTPMfile: output TPM file (by default lesionTPM.nii file in current directory)
%
% e.g. create new TPM
%    conn_createtpm('/data/109000-1_lesionmask.nii','/data/109000-1.nii');
%
% e.g. create new TPM and apply to functional dataset
%   lesionmask = '/data/109000-1_lesionmask.nii'; % lesion mask file
%   structural = '/data/109000-1.nii';            % reference structural file (in same space as lesion mask)
%   functional = '/data/109000-5.nii';            % functional data
%   newtpmfile = '/data/lesionTPM.nii';           % name of new TPM file to be created
%  
%   % creates new TPM file by adding a lesion mask to SPM's original TPM
%   conn_createtpm( lesionmask, structural, newtpmfile);
% 
%   % normalizes functional data using new TPM file
%   conn_module('preprocessing','steps','functional_segment&normalize_direct','functionals',functional,'coregtomean',false,'tpm_template',newtpmfile);
% 
%   % displays normalization results: functional + structural + gray matter mask + MNI T1 reference
%   spm_check_registration(conn_prepend('w',structural),...
%                           [conn_prepend('w',functional),',1'],...
%                           conn_prepend('wc1',structural), ...
%                           fullfile(fileparts(which('spm')),'canonical','avg152T1.nii'));


%alfnie 2020

if nargin<1||isempty(lesionmask),
    disp('Select mask file');
    [fname,fpath]=uigetfile('*.nii','Select mask file');
    if ~ischar(fpath), return; end
    lesionmask=fullfile(fname,fpath);
end
if nargin<2||isempty(structural),
    disp('Select structural file');
    [fname,fpath]=uigetfile('*.nii','Select reference structural file');
    if ~ischar(fpath), return; end
    structural=fullfile(fname,fpath);
end
if nargin<3||isempty(newTPMfile), newTPMfile=fullfile(pwd,'lesionTPM.nii'); end

% initialization
Nrepeat=2; % Number of iterations (note: set to higher value for more robust but time consuming behavior)
origTPMfile=conn_fullfile(fileparts(which('spm')),'tpm','TPM.nii');
structural=conn_file(conn_fullfile(structural)); structural=structural{1};
lesionmask=conn_file(conn_fullfile(lesionmask)); lesionmask=lesionmask{1};
newTPMfile=conn_fullfile(newTPMfile);
[nill,nill,fext]=fileparts(newTPMfile); 
if isempty(fext), newTPMfile=conn_prepend('',newTPMfile,'.nii'); end
assert(~strcmp(origTPMfile,newTPMfile),'original and target TPM files cannot be the same file'); 

% initial normalization attempt of reference structural (no masking)
fprintf('%s: initialization\n',mfilename);
conn_module('preprocessing','steps','structural_segment&normalize','structurals',structural);

Ntpm=numel(spm_vol(origTPMfile));
for nrepeat=1:Nrepeat
    fprintf('%s: iteration %d of %d\n',mfilename,nrepeat,Nrepeat);
    
    % project original mask to MNI-space
    conn_module('preprocessing','steps','structural_manualspatialdef','structurals',lesionmask,'respatialdef',conn_prepend('y_',structural));
    
    % create new TPM file with lesion mask
    newTPMfiles={};
    for n=1:Ntpm+1
        newTPMfiles{n}=conn_prepend('',newTPMfile,['.',num2str(n),'.nii']);
        if n<=Ntpm, spm_imcalc({[origTPMfile,',',num2str(n)],conn_prepend('w',lesionmask)},newTPMfiles{n},'i1.*(1-max(0,min(1,i2)))'); % other tissue classes
        else        spm_imcalc({[origTPMfile,',1'],conn_prepend('w',lesionmask)},newTPMfiles{n},'max(0,min(1,i2))');                   % lesion mask
        end
    end
    try, spm_unlink(newTPMfile); end
    spm_file_merge(newTPMfiles,newTPMfile);
    try, spm_unlink(newTPMfiles{:}); end
    
    % apply normalization to reference structural now using new TPM file
    conn_module('preprocessing','steps','structural_segment&normalize','structurals',structural,'tpm_template',newTPMfile);
end



