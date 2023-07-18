function fh=conn_image_display(varargin)
try
    fh=figure;
    if ischar(varargin{1})&&any(conn_server('util_isremotefile',varargin{1})), varargin{1}=conn_cache('pull',varargin{1}); end
    imshow(varargin{:});
    zoom on;
    set(fh,'numbertitle','off');
    if ischar(varargin{1}), set(fh,'name',varargin{1}); end
catch
    delete(fh);
    fh=[];
end
