function fh=conn_image_display(varargin)
try
    fh=figure;
    imshow(varargin{:});
    zoom on;
    set(fh,'numbertitle','off');
    if ischar(varargin{1}), set(fh,'name',varargin{1}); end
catch
    delete(fh);
    fh=[];
end
