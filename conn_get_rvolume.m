function [x,idx]=conn_get_rvolume(V,allowlinks)
persistent surf
if nargin<2||isempty(allowlinks), allowlinks=false; end

if any(conn_server('util_isremotefile',V.fname)), 
    V.fname=conn_server('util_localfile',V.fname);
    if ischar(allowlinks)&&~isempty(regexp(allowlinks,'^.*run_keepas:')), [x,idx]=conn_server('run_keepas',regexprep(allowlinks,'^.*run_keepas:',''),mfilename,V);
    elseif allowlinks, [x,idx]=conn_server('run_keep',mfilename,V);
    else [x,idx]=conn_server('run',mfilename,V);
    end
    return
end
if isempty(surf) % see CONN_gui.refs.surf
    surf.spherereduced=conn_surf_sphere(5);
    [surf.spheredefault,surf.default2reduced]=conn_surf_sphere(8,surf.spherereduced.vertices);
end

[x,idx]=conn_get_slice(V,1);
[x2,idx2]=conn_get_slice(V,conn_surf_dims(8)*[0;0;1]+1);
x=[x(:,surf.default2reduced) x2(:,surf.default2reduced)];
idx=[idx(surf.default2reduced);prod(conn_surf_dims(8))+idx2(surf.default2reduced)];



