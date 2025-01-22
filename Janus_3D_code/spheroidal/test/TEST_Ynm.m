clear;
p=8; np=(p+1)^2;
[u,v]=gl_grid(p); 
% Sns=SurfaceSph(shape_gallery(p,""));
% % Ynm test version.
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1); % scaled Gauss-legendre quadrature weights on [0,2pi](?)
wt = wt(:);

for l=1:p
    % Yl=Ynm(l,[],u,v);

    % % rescale normalization from Ynm to normal.
    % for m=-l:l
    %     y=Yl(:,m+l+1);
    %     pe=sqrt(factorial(l+abs(m))/factorial(l-abs(m))).*y;
    %     if m<0
    %         pe=(-1)^m.*pe;
    %     end
    %     Yl(:,m+l+1)=sqrt(factorial(l-m)/factorial(l+m)).*pe;
    % end

    % build from scratch using different normalization
    y = legendre2(l,cos(u))';
    yfull=[(-1).^(-l:0) ones(1,l)].*factorial(l-(-l:l))./factorial(l+(-l:l)).*y(:,[l+1:-1:2 1:l+1]);
    y =yfull.*exp(1i*v*(-l:l));
    % y = [(-1).^(-l:0) ones(1,l)].*sqrt((2*l+1)/4/pi*factorial(l-abs(-l:l))./factorial(l+abs(-l:l))).*y;
    y = sqrt((2*l+1)/4/pi*factorial(l-(-l:l))./factorial(l+(-l:l))).*y;
    Yl=y;

    err=max(max(abs(Yl'*(Yl.*repmat(wt,1,2*l+1)) - eye(2*l+1))));
    printMsg('n=%d :%4.2e\n', l, err);
end

