function pot = nonSingularDbl(S, den, opts,update_cache)
% S is a vector of vesicles of length n, DEN is a vector of D3VECs of length
% also n. OPTS is the option structure returned by SETOPTS function.
  
  persistent p pUp pDown np wt src nor iCoeff D
  
  if nargin<4,update_cache=false;end
  nv = length(S);  
  if (nv==1), pot = 0*den; return; end
  if ( isempty(wt) || p ~= opts.p)
    p     = opts.p;
    pUp   = round(opts.interUpRate * p);
    pDown = round(opts.interFilterRate * p);

    np = 2 * pUp * (pUp + 1);
    [~,gwt] = grule(pUp+1);
    wt = pi/pUp*repmat(gwt', 2*pUp, 1)./sin(parDomain(pUp));
    wt = wt(:);
  end   
  use_mat = (isempty(opts.interactUseMat) && np*nv<2e3)||opts.interactUseMat;

  toxyz = @(v) [reshape(v(1:np,:),[],1);reshape(v(np+1:2*np,:),[],1);reshape(v(2*np+1:end,:),[],1)];
  frmxyz = @(v) [reshape(v(1:nv*np,:),[],nv);reshape(v(nv*np+1:2*np*nv,:),[],nv);reshape(v(2*np*nv+1:end,:),[],nv)];
  
  %-- cache src for each time step
  if isempty(src) || update_cache
      printMsg('* Updating cached dl stokes interaction sources\n');
      src     = zeros(3*np, nv);
      nor     = zeros(3*np, nv);

      for ii=1:nv
          xUp = interpsh(S(ii).cart, pUp);
          src(:,ii) = xUp.vecForm();
          [nd3, W] = getNormal(xUp);
          nor(:,ii) = nd3.vecForm();
          iCoeff(:,ii) = (-3/4/pi)*repmat(wt .* W, 3,1);
      end
      
      if use_mat
          printMsg('* Updating cached dl stokes interaction matrix\n');
          if np*nv>3e3, warning('Constructing dl interaction matrix for %d points may take a long time and a lot of memory',np*nv);end

          %this may be buggy and/or impractical for more than 2 vesicles
          D = zeros(3*np*nv,3*np*nv);
          re = repmat(1:nv,3*np,1);re=toxyz(re);
          for iV=1:nv
              Di    = stokes_dl_smooth_matrix(toxyz(src(:,(1:nv)~=iV)),...
                  toxyz(nor(:,(1:nv)~=iV)),...
                  src(:,iV));
              rIdx  = re==iV;
              cIdx  = ~rIdx;
              D(rIdx,cIdx) = Di;
          end
      end
  end
  
  if isempty(den),return;end
    
  %-- Changing the format of data
  dblyDen = zeros(3*np, nv);
  dblyPot = zeros(3*np, nv);
  potCorr = zeros(3*np, nv);
  
  %-- Extracting points and densities from vesicles  
  for ii=1:nv    
    if(abs(S(ii).dblLayerCoeff) > 1e4*eps)
      dblyDen(:,ii) = vecForm( interpsh(den(ii), pUp) );   
      dblyDen(:,ii) = iCoeff(:,ii) .* dblyDen(:,ii);
      dblyDen(:,ii) = S(ii).dblLayerCoeff * dblyDen(:,ii);
    end
  end

  if(any(abs([S.dblLayerCoeff])>1e4*eps))
    for ii=1:nv
      if(abs(S(ii).dblLayerCoeff)>1e4*eps)
        potCorr(:,ii) = dblLayerAlltoAll(src(:,ii), dblyDen(:,ii), nor(:,ii));
      end
    end
    
    if use_mat
        dblyPot = frmxyz(D*toxyz(dblyDen));
    else
        dblyPot = dblLayerAlltoAll(src, dblyDen, nor);
        dblyPot = dblyPot - potCorr;
    end
  end
  
  dblyPot = interpsh(reshape(dblyPot,np,3*nv),p);
  dblyPot = filterSh(dblyPot, pDown);
  dblyPot = reshape(dblyPot,[],nv);

  for ii=1:nv
    pot(ii) = d3Vec(dblyPot(:,ii));
  end
