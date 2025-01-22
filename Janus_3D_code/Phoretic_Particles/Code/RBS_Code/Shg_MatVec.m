function Y = Shg_MatVec(V,DMV,L,paroffd) 

np = paroffd.np; Nb = paroffd.Nb; N = paroffd.N;   
dense=paroffd.dense; kerd=paroffd.kerd; n3=paroffd.n3; 

if strcmp(paroffd.flag_pot,'dSL_L_3D')
   mult=paroffd.eta;  
else
   mult=1; 
end

Xv = paroffd.X;   

% If V is not empty, compute A*V. Otherwise, construct A densely. 
if ~isempty(V)
    
if dense
    % Off-diagonal kernel evaluation
    K = Kernel_Eval(Xv,Xv,paroffd); 
else
    Nrv = paroffd.nor; Wv = paroffd.W2.';
    flag_pot=paroffd.flag_pot;
    pardiag=paroffd;   
    if kerd>1
       pardiag.ci = paroffd.ci(1:kerd*np); 
       pardiag.cj = paroffd.cj(1:kerd*np); 
    end
 
    YKdg = zeros(size(V)); 
    
    for k=1:n3
        idxM = (1:Nb)+Nb*(k-1); 
        pardiag.X = paroffd.X(idxM,:); 
        pardiag.nor = paroffd.nor(idxM,:); 
        pardiag.W2 = paroffd.W2(idxM); 
        YKdg(idxM,:) = Kernel_Eval(pardiag.X,pardiag.X,pardiag)*V(idxM,:); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-d + diag contributions
    if ~iscell(DMV)
       Ydg = reshape(DMV(reshape(V,Nb,[])),N,[]);
    else
      Ydg = zeros(size(V)); 
       for k=1:n3
	  ind = (1:Nb)+Nb*(k-1);
	  Ydg(ind,:) = DMV{k}(V(ind,:));  
       end  
    end
    
    if dense
    % Remove diagonal blocks
    [jj,ii] = meshgrid(1:Nb,1:Nb); 
    idb = ii+N*(jj-1); 
    id = repmat(idb(:),1,n3) + repmat(0:Nb+N*Nb:N^2,Nb^2,1); id = id(:); 
    K(id)=0; 
    
    if isempty(L)
        Y = mult*K*V + Ydg; 
    else
        Y = mult*K*V + L*V+ Ydg; 
    end
    
    else
        
        if isempty(L)
            Y = -mult*YKdg + Ydg; 
        else
            Y = -mult*YKdg + L*V+ Ydg; 
        end

	%fprintf('\n Diagonal part \n') 
	%norm(Y) 
        
        XFMM = Xv(1:kerd:end,:); 
        NrFMM = Nrv(1:kerd:end,:); 
        WFMM = Wv;
	    YFMM = RBS_FMM_Eval(full(V),WFMM,kerd,flag_pot,XFMM,XFMM,NrFMM); 
	    Y = Y + mult*YFMM; 
        %Y = Y + LOCAL_FMM_Eval(V,WFMM,kerd,flag_pot,XFMM,XFMM,NrFMM); 
        
    end
elseif dense
    % Off-diagonal kernel evaluation
    K = mult*Kernel_Eval(Xv,Xv,paroffd); 
    
    %Build dense matrix
    [jj,ii] = meshgrid(1:Nb,1:Nb); 
    idb = ii+N*(jj-1); 
    id = repmat(idb(:),1,n3) + repmat(0:Nb+N*Nb:N^2,Nb^2,1); id = id(:);
    
    if isnumeric(DMV)
        Db = repmat(DMV(:),n3,1); 
    elseif iscell(DMV)
        Db = zeros(n3*Nb^2,1); 
        for j=1:n3
           Db((1:Nb^2)+(Nb^2)*(j-1)) = DMV{j}(:);  
        end
    else
        Db = repmat(reshape(DMV(eye(Nb)),[],1),n3,1); 
    end
    
    if isempty(L)
        K(id) = Db; 
    else
        K(id) = Db+L(id); 
    end
    
    Y=K;
else
    % Return as function handle (apply)
    if ~iscell(DMV)
	DMat = @(x) DMV*x;
    else
      DMat = cell(1,n3); 
        for k=1:n3
		    DMat{k} = @(x) DMV{k}*x;
        end
    end 
    Y = @(V) Shg_MatVec(V,DMat,L,paroffd);
end

end