function [P,Y] = make_projection(p)
    %constructs projection matrix onto eigenvectors of degree 0-p;
    %NOTE: we don't project onto final column b/c of numerical reasons
    %input: p (int) degree of SpHarm
    %output: P projection matrix (p+1)^2-1 x 2*p*(p+1)
    %        Y columns of spherical harmonics
    sp =@(p) (p+1)^2;
    
    Y2=shSyn(eye(sp(p)));
    Y=Y2(:,1:sp(p)-1);
    Y(:,p^2+1)=Y(:,p^2+1)/2;
   
    [~, gwt]=g_grid(p+1);
    Sc= SurfaceSph(shape_gallery(p,''));
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);
    W = Sc.geoProp.W; W= W.*wt;
    
    W1=repmat(W,1,size(Y,2));
    
    P=Y*(W1.*Y)';
    
    %repeats matrix rep times
    %P=kron(eye(rep),P);
    
end

