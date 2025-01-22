function hOut = plotCrossSection(X, varargin)

  if(nargin<2), pauseIn = .05;end
  
  if(~isa(X,'vesicle'))
    for ii=1:length(X)
      S(ii) = vesicle(X(ii));
    end
    X = S;
  end

  for ii=1:length(X)
    C(ii,:) = X(ii).centerOfMass().vecForm()';
    XX(ii) = X(ii).cart;
  end
  
  h = plotb(XX,[],varargin{:});  
  
  C_all = mean(C,1);
  set(gca,'cameraposition', C_all);
  set(gca,'cameratarget', C_all + [0 1 0]);
  
  for ii=1:size(C,1)
    light('Position', C(ii,:));
    light('Position', C(ii,:));
    light('Position', C(ii,:));
  end
  
  if(nargout>0), hOut = h; end
