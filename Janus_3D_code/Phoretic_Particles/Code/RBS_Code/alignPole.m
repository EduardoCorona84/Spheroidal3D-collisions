function [u, v] = alignPole(T, S, pUp) 
% ALIGNPOLE() estimates the location of pole of T with respect to the
% parametrization of S, i.e. if one rotates S with the returned u and v, one
% gets to pole of the shapes aligned.
  
  if(nargin==0), testThis(); return; end

  if(nargin<3 || isempty(pUp))
    pUp = min([4*S.p 80]);
  end
  
  X = reshape(interpsh(T.cart, pUp), pUp+1, []);
  %- Finding the north pole. The coefficients can be summed at (0,0) for the
  %exact evaluation, but seems to work fine
  northPole = vec3d;
  northPole.x = mean(X.x(1,:));
  northPole.y = mean(X.y(1,:));
  northPole.z = mean(X.z(1,:));
  
  southPole = vec3d;
  southPole.x = mean(X.x(end,:));
  southPole.y = mean(X.y(end,:));
  southPole.z = mean(X.z(end,:));
  
  %- Finding the index of the closest point on the shape
  X = reshape(interpsh(S.cart, pUp), pUp+1, []);
  northPoleS = vec3d;
  northPoleS.x = mean(X.x(1,:));
  northPoleS.y = mean(X.y(1,:));
  northPoleS.z = mean(X.z(1,:));
  X = reshape(X,[],1);
            
  r = norm(X-northPole);
  rnn = norm(northPoleS-northPole);
  [~, indNorth] = min([rnn;r]);
  indNorth = indNorth-1;
  
  r = norm(X-southPole);
  rns = norm(northPoleS-southPole);
  [~, indSouth] = min([rns;r]);
  indSouth = indSouth-1;
  
  if(indNorth == 0)
    u = 0; v = 0;
    return;
  end
  
  if(indSouth == 0)
    u = pi; v = 0;
    return;
  end
  
  %- Aligning both north and south pole and check which is correct
  [u, v] = gl_grid(pUp);
  X = movePole(reshape(S.cart.to_array, [], 3), u(indNorth), v(indNorth));
  Y = movePole(reshape(S.cart.to_array, [], 3), u(indSouth), v(indSouth));
  X = vec3d(X);
  Y = vec3d(Y);
  
  errX = sum(norm(interpsh(T.cart, S.p) - X));
  errY = sum(norm(interpsh(T.cart, S.p) - Y));
  
  if(errX < errY)
    u = u(indNorth); v = v(indNorth);
  else
    u = u(indSouth); v = v(indSouth);
  end

function testThis()
  p = 8;
  shape = 'ellipseZ';
  S = SurfaceSph(shape_gallery(p, shape));
  [u, v] = gl_grid(S.p);
  thetaIdx = p/2;
  phiIdx = p/4+1;
  idx = (phiIdx-1)*(S.p+1)+thetaIdx;
  X = movePole(reshape(S.cart.to_array, [], 3), u(idx), v(idx));
  T = SurfaceSph(vec3d(X));
  figure
  subplot(2,2,1); S.plot;
  subplot(2,2,2); T.plot;
  [uR, vR] = alignPole(T, S);
  disp([uR u(idx) abs(uR-u(idx))]);
  disp([vR v(idx) abs(vR-v(idx))]);
  X = movePole(reshape(S.cart.to_array, [], 3), uR, vR);
  S.cart = vec3d(X);
  disp(max(norm(T.cart-S.cart))/max(norm(T.cart)));
  subplot(2,2,3); S.plot;
  subplot(2,2,4); T.plot;
