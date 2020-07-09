function p = conn_tcdf(x,dof)

if numel(x)==1, x=x+zeros(size(dof)); end
if numel(dof)==1, dof=dof+zeros(size(x)); end    
p=nan(size(x));
approx1=dof>1e7;
if any(approx1)
    p(approx1)=.5*erfc(-x(approx1)./sqrt(2));
end
approx2=dof<x.^2;
if any(approx2)
    p(approx2)=betainc(dof(approx2)./(dof(approx2)+x(approx2).^2),dof(approx2)/2,.5)/2;
    p(approx2&x>0)=1-p(approx2&x>0);
end
approx3=~approx1&~approx2;
if any(approx3)
    p(approx3)=.5+sign(x(approx3)).*betainc(x(approx3).^2./(dof(approx3)+x(approx3).^2),.5,dof(approx3)/2)/2;
end
