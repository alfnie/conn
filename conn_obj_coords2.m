function [oXYZ] = conn_obj_coords2(matname, XYZ)
% FORMAT [oXYZ] = conn_obj_coords2(matname,XYZ)
% Modified for SPM2 from obj_coords in SPM99. 

%Changes for SPM2  - Paul Mazaika 04/2004.
%  Reads data from sn.mat file.
%  Transformation from SPM99 -> SPM2 data formats
%       Affine -> Affine
%       Transform array reshaped by [Dims(2,:) 3] -> Tr
%          Example: 392x3 becomes 7x8x7x3, similar to Tr is 7x9x7x3.
%        **Note SPM2 uses same deformation basis functions as SPM99
%            but combined in a way to give smoother transformations.
%       [ Dims(2,:) 3 ] = size(Tr);
%       MF  ->  VF.mat
%       MG  ->  VG.mat
%       [ Dims(1,:) 2 ] -> VG.dim
%  Not used here:
%  ?    Dims(3,:)  appear to be diagonal of VF.mat. Sign change for L <-> R
%  ?    Dims(4,:)  not used
%       [ Dims(5,:) 4 ] -> VF.dim
%  ?    Dims(6,:)  appear to be diagonal of VG.mat.

if nargin < 2
  error('Not enough arguments');
end
load(deblank(matname))
tmp = find(size(XYZ)==3 | size(XYZ)==4); 
if isempty(tmp)
  % in fact can be 4 by N N by 4, to allow (fourth row = 1) format
  error('XYZ must by 3 by N, or N by 3')
elseif tmp==2
  XYZ = XYZ';
end

% from mm space to voxel space.
[x y z] = deal(XYZ(1,:),XYZ(2,:),XYZ(3,:));
Mult = inv(VG.mat);  %  was inv(MG) in SPM99
X= Mult(1,1)*x + Mult(1,2)*y + (Mult(1,3)*z + Mult(1,4));
Y= Mult(2,1)*x + Mult(2,2)*y + (Mult(2,3)*z + Mult(2,4));
Z= Mult(3,1)*x + Mult(3,2)*y + (Mult(3,3)*z + Mult(3,4));

Mult = VF.mat*Affine;   % was MF*Affine in SPM99
if (prod(size(Tr)) == 0), % no nonlinear components in sn.mat file
  X3= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*Z + Mult(1,4));
  Y3= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*Z + Mult(2,4));
  Z3= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*Z + Mult(3,4));
else % nonlinear components
  % first apply nonlinear, then affine
  %[X2,Y2,Z2] = build_transform(Transform,[Dims(2,:); Dims(1,:)],X,Y,Z);
  %   SPM2 translation:  Transform->Tr, 
  %                      Dims(2,:)->size(Tr),
  %                      Dims(1,:)->VG.dim
  %[X2,Y2,Z2] = build_transform2(Tr,[size(Tr); VG.dim],X,Y,Z);
  dim=[VG.dim,1];[X2,Y2,Z2] = build_transform2(Tr,[size(Tr); dim(1:4)],X,Y,Z);
  X3 = Mult(1,1)*X2 + Mult(1,2)*Y2 + Mult(1,3)*Z2 + Mult(1,4);
  Y3 = Mult(2,1)*X2 + Mult(2,2)*Y2 + Mult(2,3)*Z2 + Mult(2,4);
  Z3 = Mult(3,1)*X2 + Mult(3,2)*Y2 + Mult(3,3)*Z2 + Mult(3,4);

end
oXYZ = [X3;Y3;Z3];
return;

%_______________________________________________________________________
%_______________________________________________________________________
function [TX,TY,TZ] = build_transform2(T,dim,X,Y,Z)
% Note dim is size 2x4 for this SPM2 version, not 2x3 as in SPM99.
% No reshape needed. T is a 4D array.  was T = reshape(T,[dim(1,:) 3]);
TX = X;
TY = Y;
TZ = Z;
BX = basis_funk(X,dim(2,1),dim(1,1)); 
BY = basis_funk(Y,dim(2,2),dim(1,2));
BZ = basis_funk(Z,dim(2,3),dim(1,3));
for i3=1:dim(1,3),
	for i2=1:dim(1,2),
		B2 = BZ(:,:,i3).*BY(:,:,i2);
		for i1=1:dim(1,1),
			B  = B2.*BX(:,:,i1);
			TX = TX + T(i1,i2,i3,1)*B;
			TY = TY + T(i1,i2,i3,2)*B;
			TZ = TZ + T(i1,i2,i3,3)*B;
		end;
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function B = basis_funk(X,N,kk)
B = zeros([size(X) kk]);
for k=1:kk,
	if k==1,
		B(:,:,k) = ones(size(X))/sqrt(N);
	else,
		B(:,:,k) = sqrt(2/N)*cos((X-0.5)*(pi*(k-1)/N));
	end;
end;
return;
%_______________________________________________________________________
