function varargout = conn_convertcoordinates(type,varargin); %idx1,M1,dim1,M2,dim2)

switch(lower(type)),
	case 'idx2idx', % TRANSFORMS INDEXES IN ONE VOLUME TO INDEXES IN A SECOND VOLUME
		[idx1,M1,dim1,M2,dim2]=deal(varargin{:});
		
		if length(dim1)<4, dim1=[dim1,ones(1,4-length(dim1))]; end
		if length(dim2)<4, dim2=[dim2,ones(1,4-length(dim2))]; end
		[xyz1{1:4}]=ind2sub(dim1(1:4),idx1(:)); xyz1=cat(2,xyz1{:});
		xyz2=[xyz1(:,1:3),ones(size(xyz1,1),1)]*M1'*pinv(M2)';
		xyz2(:,4)=1;%xyz1(:,4);
		xyz2=max(1,min(repmat(dim2,[size(xyz2,1),1]),round(xyz2)));
		xyz2=mat2cell(xyz2,size(xyz2,1),ones(1,size(xyz2,2)));
		idx2=sub2ind(dim2,xyz2{:});
		
		varargout{1}=idx2;
		
	case 'tal2idx', % TRANSFORMS SPACE COORDINATES (mm) TO INDEXES IN A VOLUME
		[xyz1,M2,dim2]=deal(varargin{:});
		
		dim2=dim2(1:3);
		xyz2=[xyz1(:,1:3),ones(size(xyz1,1),1)]*pinv(M2)';
		xyz2=xyz2(:,1:3);
		xyz2=max(1,min(repmat(dim2,[size(xyz2,1),1]),round(xyz2)));
		xyz2=mat2cell(xyz2,size(xyz2,1),ones(1,size(xyz2,2)));
		idx2=sub2ind(dim2,xyz2{:});
		
		varargout{1}=idx2;
		
	case 'idx2tal', % TRANSFORMS INDEXES IN A VOLUME TO SPACE COORDINATES (mm)
		[idx1,M1,dim1]=deal(varargin{:});
		
		[xyz1{1:4}]=ind2sub(dim1(1:3),idx1(:)); xyz1=cat(2,xyz1{:});
		xyz1=[xyz1(:,1:3),ones(size(xyz1,1),1)]*M1';
		
		varargout{1}=xyz1;
		
end
	