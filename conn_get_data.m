function [X,IDX,Q]=conn_get_data(V,xyz,weight,components)
% CONN_GET_DATA Reads data from file
%    [X,IDX]=CONN_GET_DATA(V,XYZ) reads data from volume V at coordinates XYZ 
%	   V is a volume structure returned by CONN_VOL
%	   XYZ is a Nx3 matrix of voxel positions (in mm)
%	   X is a data matrix of size MxK (time-points by voxels), where K<=N (repeated and missing voxels in the volume space are discarded)
%	   and IDX is an index vector of size Nx1 such that if idx(n)>0 then X(:,idx(n)) is the time-series at position XYZ(n,:) 
%	   -note: IDX is zero for missing data
%    [X,IDX]=CONN_GET_DATA(V,XYZ,SAMPLE) reads one time-sample data from volume V at coordinates XYZ at sample number SAMPLE 
%
%    X=CONN_GET_DATA(V,XYZ,WEIGHT,COMPONENTS) obtains a number of characteristic time-series
%		from the data (average time series, if COMPONENTS=1; or average+PCA components,
%		if COMPONENTS>1) where each voxel activation is optionally weighted by the values
%		in vector WEIGHT (a vector of N values).
%	   X is a data matrix of size MxCOMPONENTS (time-points by number of components) 
%


if nargin<3, weight=[]; end
if nargin<4 || isempty(components), if length(weight)<=1, components=inf; else, components=1; end; end
if any(conn_server('util_isremotefile',V.fname)), V.fname=conn_server('util_localfile',V.fname); [X,IDX,Q]=conn_server('run',mfilename,V,xyz,weight,components); return; end

if ~components, X=[]; return; end
if isfield(V,'softlink')&&~isempty(V.softlink), 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    matcfilename=fullfile(file_path,V.softlink); 
else
    matcfilename=[V.fname,'c'];
end
    
idx=conn_convertcoordinates('tal2idx',xyz,V.matdim.mat,V.matdim.dim);
[idx,idxbak,idx0]=unique(idx);				
[nill,voxels,idx1]=intersect(V.voxels,idx);
[voxels,idx2]=sort(voxels); 
temp=zeros(1,length(idx)); temp(idx1(idx2))=(1:length(idx2)); IDX=temp(idx0); 
idxbak=idxbak(idx1(idx2));

if isinf(components), % extracts activity timecourse for each of the voxels of the ROI
	if isempty(weight), X=nan(V.size.Nt,length(voxels));
	else, X=zeros(1,length(voxels)); end
	n0=1;
	handle=fopen(matcfilename,'rb');
	for n1=1:length(voxels),
		if voxels(n1)>n0,fseek(handle,4*(voxels(n1)-n0)*V.size.Nt,0);end
		if isempty(weight), X(:,n1)=fread(handle,V.size.Nt,'float');
		else, fseek(handle,4*(weight-1),0); X(n1)=fread(handle,1,'float',4*(V.size.Nt-weight)); end
		n0=voxels(n1)+1;
	end
	fclose(handle);
	X=V.scale*X;
elseif components==1, % extracts (weighted) average ROI activation
	if isempty(weight), weight=ones(size(voxels))/length(voxels);
	else, weight=weight(idxbak)/max(eps,sum(weight(idxbak))); end
	X=zeros(V.size.Nt,1);
	n0=1;
	handle=fopen(matcfilename,'rb');
	for n1=1:length(voxels),
		if voxels(n1)>n0,fseek(handle,4*(voxels(n1)-n0)*V.size.Nt,0);end
		X=X+fread(handle,V.size.Nt,'float')*weight(n1);
		n0=voxels(n1)+1;
	end
	fclose(handle);
	X=V.scale*X;
else, % returns (weighted) average+PCA ROI activation
	if isempty(weight), weight=ones(size(voxels))/length(voxels);
	else, weight=weight(idxbak)/max(eps,sum(weight(idxbak))); end
	M=zeros(V.size.Nt,1);
	X=zeros(V.size.Nt,V.size.Nt);
	n0=1;
	handle=fopen(matcfilename,'rb');
	for n1=1:length(voxels),
		if voxels(n1)>n0,fseek(handle,4*(voxels(n1)-n0)*V.size.Nt,0);end
		x=fread(handle,V.size.Nt,'float');
		M=M+weight(n1)*x;
		X=X+weight(n1)*x*x';
		n0=voxels(n1)+1;
	end
	fclose(handle);
	M=V.scale*M;
	X=V.scale*V.scale*X;
	m=M/max(eps,norm(M)); 
	X=X-X*m*m'-m*m'*X+m*m'*X*m*m';
	[Q,D]=svd(X);
	X=[M,Q(:,1:components-1)*diag(sqrt(abs(diag(D(1:components-1,1:components-1)))))];
end
