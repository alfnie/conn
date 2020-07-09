function conn_premakestandalone(varargin)
% prepares CONN for standalone compilation
% see conn/standalone/conn_makestandalone to compile CONN
% or run conn_premakestandalone before creating your own SPM+toolboxes compilation

for f={'conn_batch','conn_grid'},
    fname=which(f{1});
    str=help(fname);
    fh=fopen(regexprep(fname,'\.m$','.help.txt'),'wt');
    fprintf(fh,'\n%s\n',str);
    fclose(fh);
end
