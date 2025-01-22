function cent = getCenter(S)

persistent wt p

if(isempty(wt) || p ~= S(1).p)
  p = S(1).p;
  [~, gwt]=g_grid(p+1);
  wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
  wt = wt(:)';
end

for ii=1:length(S)
  W = S(ii).geoProp.W;
  X = S(ii).cart;
  N = S(ii).geoProp.nor;

  V = 1/3*wt*(dot(X,N).*W);

  cx(ii,1) = wt * (X.x(:).^2 .* N.x .*W)/2/V;
  cy(ii,1) = wt * (X.y(:).^2 .* N.y .*W)/2/V;
  cz(ii,1) = wt * (X.z(:).^2 .* N.z .*W)/2/V;
end
cent = vec3d([cx cy cz]);
  
