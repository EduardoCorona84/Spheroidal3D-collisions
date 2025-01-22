function I = integrateOverS(S, f)

  persistent wt p
  
  if(isempty(wt) || p ~= S.p)
    p = S.p;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
  end

  W = S.geoProp.W;
  if(isa(f,'vec3d'))
    I = vec3d;
    I.x = wt*(W.*f.x);
    I.y = wt*(W.*f.y);
    I.z = wt*(W.*f.z);
  else
    I = wt*(W.*f);
  end
