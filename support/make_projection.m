function [P,Y] = make_projection(p,rep)
    %constructs projection matrix onto eigenvectors of degree 0-p;
    %NOTE: we don't project onto final column b/c of numerical reasons
    %input: p (int) degree of SpHarm
    %output: P projection matrix (p+1)^2-1 x 2*p*(p+1)
    %        Y columns of spherical harmonics
    if nargin <2
        rep=1;
    end
    np = @(p) 2*(p+1)*p;
    sp =@(p) (p+1)^2;
    [u v]=gl_grid(p);
    ii=(1:sp(p))';
    nn=floor(sqrt(ii-1));
    mm=ii-nn.^2-nn-1;
   
    Y2=shSyn(eye((p+1)^2));
    Y=Y2(:,1:(p+1)^2-1);
    Y(:,p^2+1)=Y(:,p^2+1)/2;
  
   
    [~, gwt]=g_grid(p+1);
    Sc= SurfaceSph(shape_gallery(p,''));
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);
    W = Sc.geoProp.W; W= W.*wt;
    
  
    W1=repmat(W,1,size(Y,2));
    W2=repmat(W,1,((2*p*(p+1))));
    size(W2)
    
    P=Y*(W1.*Y)';
    D=diag((W1.*Y)'*Y);

    %repeats matrix rep times
    P=kron(eye(rep),P);
    
end

