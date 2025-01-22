function [P,V] = make_vshprojection(p,perm)
    if nargin==1
        perm=false; 
    end

    %input: p (int) degree of SpHarm
    %output: P projection matrix 6*p*(p+1) x 6*p*(p+1)
    %V columns of vector spherical harmonics
    sp =@(p) (p+1)^2;
    ii = 1:sp(p); 
    nn = floor(sqrt(ii-1)); 
    nrm = [(2*nn+1).*(nn+1) (2*nn+1).*nn nn.*(nn+1)];
    
    V=VshSyn(eye(3*sp(p)),'VW');
    dl = [(sp(p-1)+1):(sp(p)+1) (sp(p)+sp(p-1)+1):(2*sp(p)+1) (2*sp(p)+sp(p-1)+1):(3*sp(p))];
    V(:,dl)=[]; 
    nrm(dl)=[];
    
    [~, gwt]=g_grid(p+1);
    Sc= SurfaceSph(shape_gallery(p,''));
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);
    W = Sc.geoProp.W; W= W.*wt;
    
    W1=repmat(W,3,size(V,2)); 
    N1 = repmat(1./nrm,size(V,1),1); 
    
    P=V*(W1.*V.*N1)';
    
    % permute rows of projection matrix
    if perm
       np = 2*p*(p+1); 
       Jprm = [1:3:3*np 2:3:3*np 3:3:3*np]; 
       P(Jprm,Jprm) = P; 
    end
    
end