function [prec,prev] = RBS_setprec(ID,n3,prtype,p,params,prev)

if strcmp(prtype,'bkdiag') || nargin<4
    if ~iscell(ID)
       Nb = size(ID,2); 
       N = Nb*n3;    
       prec = @(V) reshape(ID*reshape(V,Nb,[]),N,[]);  
    else
       Nb = size(ID{1},2); 
       N = n3*Nb; 
       JT=zeros(n3*Nb^2,1); 
        IT=JT; VT=JT; 
        for k=1:n3
            idxM = (1:Nb)+Nb*(k-1);
            idJ = (1:Nb*Nb)+Nb*Nb*(k-1); 
            [JJ,II] = meshgrid(idxM,idxM); 
            %VV = ITSSDd{k}(:);
            JT(idJ)=JJ(:); 
            IT(idJ)=II(:); 
            VT(idJ)=ID{k}(:);  
        end
        SITDd = sparse(IT,JT,VT,N,N); 
        prec = @(V) SITDd*V; 
    end
    
    prev=[]; 
    
elseif strcmp(prtype,'TT')
    if isempty(prev)
        prev=cell(2,1); 
    end
    
    params.NP = [3 factor(p+1) factor(p) 2];  
    Kfun = @(ind) LOCAL_Tensor_Kernel_Eval(ind,ID,params,'dense'); 
    fprintf('\n TT preconditioner build \n')
    % Tensor dimensions
    Lg = log2(n3); 
    row = [params.NP 2*ones(1,Lg)]; col=row; 
    RC = row.*col; 
        
    % AMEN parameters
    acc = 1e-2; 
    params.y0 = prev{1}; 
    params.acc = 0.1*acc;        
    params.rmax = 200; 
    if isempty(prev{1})
        params.kr = 32;
    else
        params.kr=16; 
    end
    params.nswp = ceil(params.rmax/params.kr)+2;       
    params.exit_acc = acc;   
        
    % TT compression and inversion
    fprintf('\n TT compression \n')  
    TTK = amen_cross(RC,Kfun,params.acc,'verb',1,'vec',1,'nswp',params.nswp,...
        'kickrank',params.kr,'tol_exit',params.exit_acc,'y0',params.y0);        
    TTK = tt_matrix(round(TTK,acc),row,col);    
    display(TTK)
        
    fprintf('\n TT inversion \n')             
    params.y0 = prev{2}; 
    if isempty(prev{1})
        params.kr = 32;
    else
        params.kr=8; 
    end
    params.ls = 'gmres'; params.pr = [];  
    params.nmax=2000;             
        
    %TTkron = kron2(TTK,tt_eye(TTK.n));  
    %TTKM = TTkron;
    TTI = get_tt_inverse(TTK,params);      
    display(TTI)   
        
    prec = @(x) mtimes_ec(TTI,x);   
    prev{2} = tt_tensor(TTI); 
    prev{1} = tt_tensor(TTK); 
    
end

end

function A = LOCAL_Tensor_Kernel_Eval(ind,D,par,type)


nlp = length(par.NP); 

shIdx   = ind(:,1:nlp)-1;
centIdx = ind(:,nlp+1:end)-1;

if size(ind,2)>nlp 
    I = mod(centIdx,2);
    J = (centIdx-I)/2;    
    I = bi2de(I);
    J = bi2de(J);
else
    I = []; J = []; 
end

if nlp>1
    rNP = repmat(par.NP,size(ind,1),1);
    pNP = prod(par.NP); 
    
    IP = mod(shIdx, rNP); 
    JP = (shIdx-IP)./rNP;
    
    pw = repmat(cumprod([1 par.NP(1:end-1)]),size(ind,1),1);
    IP = sum((pw.*IP).').';
    JP = sum((pw.*JP).').';     
else
    pNP = par.NP; 
    IP = mod(shIdx, par.NP); 
    JP = (shIdx-IP)./par.NP;
end

Dist = abs(I-J);
isdiag = ~(Dist>0);     
   
I = I*pNP;        
J = J*pNP;  

if ~isempty(I) 
   II = IP+I+1; 
   JJ = JP+J+1;    
else
   II = IP+1; 
   JJ = JP+1;    
end

A = zeros(size(ind,1),1); 

sdg = sum(isdiag); 

% Diagonal block entries (may change to a sparse matrix format)
if strcmp(type,'dense')
    if ~iscell(D)
        A(isdiag) = D(IP(isdiag)+1+pNP*JP(isdiag));
    elseif sdg>0
        idx=1:size(ind,1); 
        idd = idx(isdiag); 
        for k=1:sdg
           A(idd(k)) = D{I(idd(k))/pNP+1}(IP(idd(k))+1+pNP*JP(idd(k)));        
        end
    end
elseif sdg>0
    D = par.S; 
    cache = (5000)^2;  
    
    % Rotfree (use LkernelD_rows)
    IId = II(isdiag); JJd = JJ(isdiag); 
    
    N = pNP*(2^size(centIdx,2));  
    Atmp = D(IId+N*(JJd-1)); 
    ind0 = Atmp==0; 
    nzr = sum(ind0);
    nnzr = sum(~ind0); 
    [x,~,~] = find(D);   
     
    fprintf('S nnz %d:%d,%d \n',length(x),full(nnzr),full(nzr));   
    
    if nzr>0
    % LOCAL rows and columns
    IPd = IP(isdiag); JPd = JP(isdiag);
    IPd = IPd(ind0); JPd = JPd(ind0); 
    % Block index
    Id = I(isdiag); 
    Bd = Id./pNP+1; Bd = Bd(ind0); 
    
    [Bs,Js] = sort(Bd); 
    IPs = IPd(Js); 
    
    cuts = 1:sum(ind0);  
    
    if nzr>1
        cuts = cuts((Bs(2:end)-Bs(1:end-1))>0);
    else
        cuts = zeros(1,0);   
    end
    
    nsurf = length(cuts)+1; 
    cuts = [0 cuts nzr]; 
    Surf = cell(nsurf,1); rows = Surf; 
    
    for k=1:nsurf
       rows{k} = IPs(cuts(k)+1:cuts(k+1))+1;  
       Surf{k} = par.Surf{Bs(cuts(k)+1)}; 
    end
    
    % Find rows of blocks using LkernelD_rows
    Arows = LkernelD_rows(rows,Surf,Js,'rotfree');  
    % Find the values we're looking for
    indtmp = (JPd.'+1 + pNP*(0:nzr-1));  
    a = -0.5; Diag = (repmat(IPd+1,1,pNP)==repmat(1:pNP,nzr,1));       
    Atmp(ind0) = a*Diag(indtmp)+Arows(indtmp)+(0.5/pNP)*ones(1,nzr);    
    
    % Incorporate rows to D.  
    [Di,Dj,Dv] = find(D); 
    % Rows and Columns of original matrix
    [jj,ii] = meshgrid(1:pNP,IId(ind0)); 
    jj = jj + repmat(Id(ind0),1,pNP); 
    % Add to existing rows, cols and values
    Di = [Di ; ii(:)]; Dj = [Dj ; jj(:)]; 
    Dv = [Dv ; a*Diag(:)+Arows(:)+(0.5/pNP)*ones(size(ii(:)))];
    % Discard repeated entries (sparse adds repeats) 
    U = unique([Di Dj Dv],'rows','first'); 
    
    if size(U,1)>cache   
       cnr = length(unique(ii(:)));        
       U = U(find(~ismember(U(:,1),par.rowlist(1:cnr))),:);     
       par.rowlist(1:cnr) = [];  
    end
    
    % Build new sparse matrix
    D = sparse(U(:,1),U(:,2),U(:,3),N,N);       
    
    if strcmp(parcase,'params')
        params.S = D; 
        params.rowlist = [par.rowlist ; unique(ii(:))];     
    elseif strcmp(parcase,'parblk')
        parblk.S = D;  
        parblk.rowlist = [par.rowlist ; unique(ii(:))]; 
    end 
    end
    
    A(isdiag) = Atmp;   
end

% Off-diagonal block entries
isoffd = ~isdiag; 
if sum(isoffd)>0
par.W2 = par.W2(JJ(isoffd)); 
par.nor = par.nor(JJ(isoffd),:);    
X = par.X(II(isoffd),:);
Y = par.X(JJ(isoffd),:);
if isfield(par,'ci')
   par.ci = par.ci(II(isoffd)); 
   par.cj=par.cj(JJ(isoffd));  
end

Aoff = Kernel_Eval_vec(X,Y,par);   
A(isoffd) = Aoff;
end     

end