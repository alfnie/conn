function h=conn_waitbar(varargin)
global CONN_x;
h=[];
if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display),
    if ischar(varargin{1}), 
        if strcmp(varargin{1},'redraw'), h=varargin{2}; if ishandle(h), figure(h); end
        else h=varargin{2}; delete(h(ishandle(h)));
        end
    else h=conn_timedwaitbar(varargin{:}); end
else
    if ischar(varargin{1}), 
        if strcmp(varargin{1},'redraw'),
        else conn_cumdisp;
        end
    elseif isnan(varargin{1}), 
    elseif ~varargin{1}, 
        conn_cumdisp; conn_disp(varargin{2}); 
    else
        str=[num2str(100*varargin{1},'%3.1f'),'% '];
        if nargin>2,str=[str, ' (',varargin{3},')']; end
        conn_cumdisp(str); 
    end
end
end
function conn_cumdisp(txt)
persistent oldtxt;

if nargin<1,
    oldtxt=''; 
    fprintf(1,'\n'); 
else
    fprintf(1,[repmat('\b',[1,length(oldtxt)]),'%s'],txt);
    oldtxt=sprintf('%s',txt);
end
end
