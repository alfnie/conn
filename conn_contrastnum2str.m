function str = conn_contrastnum2str(num)

try
    str=regexprep(cellstr(rats(num)),'\s+',' ');
    str=regexprep(sprintf('%s;',str{:}),';$','');
    tol = max(abs(reshape(str2num(str),1,[]) - reshape(num,1,[])));
    if tol>1e-6, str=mat2str(num); end
catch
   str=mat2str(num);
end

    