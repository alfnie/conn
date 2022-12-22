function [q,p] = conn_signflip(files,w)
% internal
% flips signs of vector field in the direction of maximal variance

if nargin<2, w=[]; end
if isnumeric(files), % conn_singflip(x) - enter raw data directly (observations-by-samples matrix [Nd x N1 x N2 x ...])
    q=files; 
    mask=shiftdim(any(q~=0,1),1);
    if ~isempty(w), mq=w;
    else
        mq0=mean(q(:,mask),2);
        [mq,d]=svd(q(:,mask)*q(:,mask)');
        mq=mq(:,1);
        if mq'*mq0<0, mq=-mq; end % first singular vector (with sign/direction fixed to that of the mean vector) 
    end
    p=zeros(size(mask));
    p(mask)=mq'*q(:,mask);
    q(:,p<0)=-q(:,p<0);
else % conn_signflip(files,w)
    if any(conn_server('util_isremotefile',files)), [varargout{1:nargout}]=conn_server('run',mfilename,conn_server('util_localfile',files),w); return; end    
    files=conn_server('util_localfile',files);
    if iscell(files), files=char(files); end
    if ~isstruct(files), files=spm_vol(files); end 
    p=0;
    for n=1:numel(files)
        q=spm_read_vols(files(n));
        p=p+w(n)*q;
    end
    doflip=p<0;
    if nnz(doflip)
        for n=1:numel(files)
            q=spm_read_vols(files(n));
            q(doflip)=-q(doflip);
            spm_write_vol(files(n),q);
        end
    end
end


