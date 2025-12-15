function Y = Vsh_Kernel_Eval_off(F,type,doAna,out,rho,u,v,Nr)
%{
This code evaluates the integral kernel given by type and density given by
F on spherical coordinates (rho,u,v). 

Inputs
F - (double 3*(np / sp) x nv) integral density or corresponding vector
spherical harmonics coefficients 
type - kernel to be evaluated (SMat, SpMat, DMat, TMat, etc). 
doAna - (bool) indicates whether F is density or spherical harmonic coefs
out - (bool) exterior vs interior problem
(rho,u,v) - (double ntrg x 1) spherical coordinates arrays. Code evaluates
for rho=1, [u v] = gl_grid(p) if empty using fast transforms. 
Nr - (double ntrg x 3) Normal vector (if needed). 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent ShTg_mem ShTr_mem ShY_mem Sc_mem
if isempty(ShTg_mem)
    ShTg_mem = containers.Map('KeyType', 'double', 'ValueType', 'any');
    ShTr_mem = containers.Map('KeyType', 'double', 'ValueType', 'any');
    ShY_mem = containers.Map('KeyType', 'double', 'ValueType', 'any');
    Sc_mem = containers.Map('KeyType', 'double', 'ValueType', 'any');
end
if nargin ==1 
    % in this case F is p and we are just initializng memory
    p = F;
    if ~ShTg_mem.isKey(p)
        [ShTg_mem(p),ShTr_mem(p),ShY_mem(p)] = get_Surfgrad_coeffs(ceil(1.5*p),p);
        Sc_mem(p) = SurfaceSph(shape_gallery(p,''));
    end
    return 
end
%% F is 3*(np/sp) x nvecs
[d1,d2] = size(F);
%% Accept F or Fh as inputs
if d1>1
    if doAna
        np = d1/3;
        p = (sqrt(2*np+1)-1)/2;
        Vshc = VshAna(F,'VW');
    else
        sp = d1/3;
        p = sqrt(sp)-1;
        Vshc = F;
    end
else
    p = F(1); d2 = F(2);
    if size(F,2)>2
        Rdiag = F(3:end).';
    else
        Rdiag=[];
    end
end
%% Coefficient matrices (add option to load only once)
try 
    ShTg = ShTg_mem(p); ShTr = ShTr_mem(p); ShY = ShY_mem(p); Sc = Sc_mem(p);
catch
    [ShTg_mem(p),ShTr_mem(p),ShY_mem(p)] = get_Surfgrad_coeffs(ceil(1.5*p),p);
    Sc_mem(p) = SurfaceSph(shape_gallery(p,''));
    ShTg = ShTg_mem(p); ShTr = ShTr_mem(p); ShY = ShY_mem(p); Sc = Sc_mem(p);
end
if d1>1
    % Synthesize
    if d2>1 || ~isempty(u)
        %% LOCAL_off_VshSyn
        [~,d2] = size(Vshc);
        ntrg = 3*length(u)/d2;
        %np = 2*p*(p+1);
        sp = (p+1)^2;
        ii = (1:sp)';
        nn = floor(sqrt(ii-1));
        Y = zeros(ntrg,1);
        er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];

        %R = repmat(rho,3,sp);
        %R = repmat(reshape(repmat(rho.',3,1),[],1),1,sp);
        %P = repmat(nn.',ntrg*d2,1);
        R = repmat(rho,1,p+1);
        P = repmat(0:p,length(rho),1);
        iir = reshape(repmat((1:length(rho)),3,1),[],1);
        iip = nn+1;
        rpf = @(x,i,j) x(i,j);

        if ~strcmp(type(1),'T')
            % Get functions for coefficients for V,W and X functions
            [FV,FW,FX] = LOCAL_get_Frn(type,out);
        else
            % Get functions for coefficients for V,W and X functions
            lay =  [type(2) 'Mat']; dlay = [type(2) 'pMat'];
            [FV,FW,FX] = LOCAL_get_Frn(dlay,out);
        end

        if out
            %Exterior
            Fr = [rpf(FV(R,P),iir,iip) rpf(FW{2}(R,P),iir,iip) rpf(FW{1}(R,P),iir,iip) rpf(FX(R,P),iir,iip)];
            qhex = [Vshc(1:sp,:);Vshc(sp+1:2*sp,:);Vshc(sp+1:end,:)];
        else
            %Interior
            Fr = [rpf(FV{1}(R,P),iir,iip) rpf(FV{2}(R,P),iir,iip) rpf(FW(R,P),iir,iip) rpf(FX(R,P),iir,iip)];
            qhex = [Vshc(1:sp,:);Vshc(1:sp,:);Vshc(sp+1:end,:)];
        end

        if ~strcmp(type(1),'T')
            %Sqh = [SVext(R,nn) ; SWVext(R,nn) ; SWWext(R,nn) ; SXext(R,nn)].*qhex;
            Mqh = Fr.*(reshape(repmat(qhex,ntrg,1),4*sp,[]).');
        else
            % Get functions for Single / Double layer for V,W and X functions
            [SV,SW,SX] = LOCAL_get_Frn(lay,out);
            if out
                %Layer potential f(r,n)/r
                SFr = [rpf(SV(R,P)./(R+(R==0)),iir,iip) rpf(SW{2}(R,P)./(R+(R==0)),iir,iip) ...
                    rpf(SW{1}(R,P)./(R+(R==0)),iir,iip) rpf(SX(R,P)./(R+(R==0)),iir,iip)];
                %Pressure
                PrFr = rpf(FW{3}(R,P),iir,iip);
                MPr = PrFr.*(reshape(repmat(Vshc(sp+1:2*sp,:),ntrg,1),sp,[]).');

                %indices
                idpr = sp+1:2*sp;
                idx = [1:sp ; 1:sp ; (sp+1):2*sp ; (2*sp+1):3*sp]; idx=idx';
            else
                %Layer potential f(r,n)/r
                SFr = [rpf(SV{1}(R,P)./(R+(R==0)),iir,iip) rpf(SV{2}(R,P)./(R+(R==0)),iir,iip) ...
                    rpf(SW(R,P)./(R+(R==0)),iir,iip) rpf(SX(R,P)./(R+(R==0)),iir,iip)];
                %Pressure
                PrFr = rpf(FV{3}(R,P),iir,iip);
                MPr = PrFr.*(reshape(repmat(Vshc(1:sp,:),ntrg,1),sp,[]).');

                %indices
                %idpr = 1:sp;
                idpr = sp+1:2*sp;
                idx = [1:sp ; (sp+1):2*sp ; (sp+1):2*sp ; (2*sp+1):3*sp]; idx=idx';
            end

            eu = repmat(sin(u),1,3).*[cos(u).*cos(v) cos(u).*sin(v) -sin(u)];
            ev = repmat(sin(u),1,3).*[-sin(v) cos(v) zeros(size(u,1),1)];

            Ndotr = sum(Nr(1:3:end,:).*er,2);
            Ndotr = repmat(Ndotr.',3,1);
            Ndotu = sum(Nr(1:3:end,:).*eu,2)./sin(u).^2;
            Ndotu = repmat(Ndotu.',3,1);
            Ndotv = sum(Nr(1:3:end,:).*ev,2)./sin(u).^2;
            Ndotv = repmat(Ndotv.',3,1);

            Ndote = [Ndotu(:) Ndotv(:) Ndotr(:)];

            % Radial and Tangential contributions
            MRqh = Fr.*qhex.';
            MGqh = SFr.*qhex.';
            
            Mqh = zeros(ntrg*d2,3*sp);
            for j=1:3
                Ndote_j = Ndote(:,j);
                % These are sparse so transpose is fast 
                % turns out concatenating the sparse matrices is worth the memory overhead
                ShTg_j_T = cat(2, ShTg{j}(:,idx(:,1)), ShTg{j}(:,idx(:,2)), ShTg{j}(:,idx(:,3)), ShTg{j}(:,idx(:,4))).';
                ShTr_j_T = cat(2, ShTr{j}(:,idx(:,1)), ShTr{j}(:,idx(:,2)), ShTr{j}(:,idx(:,3)), ShTr{j}(:,idx(:,4))).';
                ShY_ji_T = ShY{j}(:,idpr).';
                % Now do the mat-mat all at once for speed of computation
                tmp1 = MGqh*ShTg_j_T;
                tmp2 = MRqh*ShTr_j_T;
                tmp3 = MPr*ShY_ji_T;
                Mqh = Mqh + Ndote_j.*(tmp1+tmp2+tmp3);
            end
        end

        Vk = Vnm('VWX',Sc,0:p,[],u,v,er);

        if size(Mqh,2)==3*sp
            for k=0:p
                ind = k^2+1:(k+1)^2;
                indv = [ind ind+sp ind+2*sp];
                Yk = sum([Vk{1}(:,ind) Vk{2}(:,ind) Vk{3}(:,ind)].*Mqh(:,indv),2);
                Y = Y + sum(reshape(Yk,[],d2),2);
            end
        else
            for k=0:p
                ind = k^2+1:(k+1)^2;
                indv = [ind ind+sp ind+2*sp ind+3*sp];

                if out
                    %Yk = sum([Vk Vk Wk Xk].*Mqh(:,indv),2);
                    Yk = sum([Vk{1}(:,ind) Vk{1}(:,ind) Vk{2}(:,ind) Vk{3}(:,ind)].*Mqh(:,indv),2);
                else
                    %Yk = sum([Vk Wk Wk Xk].*Mqh(:,indv),2);
                    Yk = sum([Vk{1}(:,ind) Vk{2}(:,ind) Vk{2}(:,ind) Vk{3}(:,ind)].*Mqh(:,indv),2);
                end

                Y = Y + sum(reshape(Yk,[],d2),2);
            end
        end
        return 
    else
        %% LOCAL_self_VshSyn
        sp = (p+1)^2;
        ii = (1:sp)';
        nn = floor(sqrt(ii-1));
        R = 1; P = nn;

        if ~strcmp(type(1),'T')
            % Get functions for coefficients for V,W and X functions
            [FV,FW,FX] = LOCAL_get_Frn(type,out);
        else
            % Get functions for coefficients for V,W and X functions
            lay =  [type(2) 'Mat']; dlay = [type(2) 'pMat'];
            [FV,FW,FX] = LOCAL_get_Frn(dlay,out);
        end

        if out
            %Exterior
            Fr = [FV(R,P);FW{2}(R,P);FW{1}(R,P);FX(R,P)];
            qhex = [Vshc(1:sp);Vshc(sp+1:2*sp);Vshc(sp+1:end)];
        else
            %Interior
            Fr = [FV{1}(R,P);FV{2}(R,P);FW(R,P);FX(R,P)];
            qhex = [Vshc(1:sp);Vshc(1:sp);Vshc(sp+1:end)];
        end

        if ~strcmp(type(1),'T')
            Mqh = Fr.*qhex;

            if out
                Mqh = [Mqh(1:sp)+Mqh(sp+1:2*sp) ; Mqh(2*sp+1:end)];
            else
                Mqh = [Mqh(1:sp); Mqh(2*sp+1:3*sp)+Mqh(sp+1:2*sp) ; Mqh(3*sp+1:end)];
            end
        else
            % Get functions for Single / Double layer for V,W and X functions
            [SV,SW,SX] = LOCAL_get_Frn(lay,out);
            if out
                %Layer potential f(r,n)/r
                SFr = (1/R)*[SV(R,P);SW{2}(R,P);SW{1}(R,P);SX(R,P)];
                %Pressure
                PrFr = FW{3}(R,P);
                MPr = PrFr.*Vshc(sp+1:2*sp);

                %indices
                idpr = sp+1:2*sp;
                idx = [1:sp ; 1:sp ; (sp+1):2*sp ; (2*sp+1):3*sp]; idx=idx';
            else
                %Layer potential f(r,n)/r
                SFr = (1/R)*[SV{1}(R,P);SV{2}(R,P);SW(R,P);SX(R,P)];
                %Pressure
                PrFr = FV{3}(R,P);
                MPr = PrFr.*Vshc(1:sp);

                %indices
                idpr = sp+1:2*sp;
                %idpr = 1:sp;
                idx = [1:sp ; (sp+1):2*sp ; (sp+1):2*sp ; (2*sp+1):3*sp]; idx=idx';
            end

            % Radial and Tangential contributions
            MRqh = Fr.*qhex;
            MGqh = SFr.*qhex;
            Mqh = zeros(3*sp,1);

            for i=1:4
                idq = (1:sp) + (i-1)*sp;
                Mqh = Mqh + ShTg{3}(:,idx(:,i))*MGqh(idq) + ShTr{3}(:,idx(:,i))*MRqh(idq);
            end

            Pqh = ShY{3}(:,idpr)*MPr;
            Mqh = Mqh + Pqh;
        end

        Y = VshSyn(Mqh,'VW');
        Y = reshape(reshape(Y,[],3).',[],1);
        return
    end
else
    %% Build dense matrix
    if isempty(Sc)
        Sc = SurfaceSph(shape_gallery(p,''));
    end

    ntrg = 3*length(u)/d2;
    np = 2*p*(p+1);
    sp = (p+1)^2;
    ii = (1:sp)';
    nn = floor(sqrt(ii-1));
    er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];

    Jprm = [1:3:3*np 2:3:3*np 3:3:3*np];
    Ji(Jprm) = 1:3*np;
    Idt = eye(3*np);
    T = VshAna(Idt(:,Ji),'VW'); %Analysis operator

    R = repmat(reshape(repmat(rho.',3,1),[],1),1,sp);
    P = repmat(nn.',ntrg*d2,1);

    if ~strcmp(type(1),'T')
        % Get functions for coefficients for V,W and X functions
        [FV,FW,FX] = LOCAL_get_Frn(type,out);
    else
        % Get functions for coefficients for V,W and X functions
        lay =  [type(2) 'Mat']; dlay = [type(2) 'pMat'];
        [FV,FW,FX] = LOCAL_get_Frn(dlay,out);
    end

    if out
        %Exterior
        Fr = [FV(R,P) FW{2}(R,P) FW{1}(R,P) FX(R,P)];
        idqh = [ (1:sp).' ; (sp+1:2*sp).' ; (sp+1:3*sp).'];
    else
        %Interior
        Fr = [FV{1}(R,P) FV{2}(R,P) FW(R,P) FX(R,P)];
        idqh = [ (1:sp).' ; (1:sp).' ; (sp+1:3*sp).'];
    end

    if ~strcmp(type(1),'T')
        Mqh = Fr;
    else
        % Get functions for Single / Double layer for V,W and X functions
        [SV,SW,SX] = LOCAL_get_Frn(lay,out);
        if out
            %Layer potential f(r,n)/r
            SFr = [SV(R,P) SW{2}(R,P) SW{1}(R,P) SX(R,P)]./repmat(R+(R==0),1,4);
            %Pressure
            PrFr = FW{3}(R,P);
            MPr = PrFr; %.*(reshape(repmat(Vshc(sp+1:2*sp,:),ntrg,1),sp,[]).');

            %indices
            idprin = sp+1:2*sp;
            idpr = idprin;
            idx = [1:sp ; 1:sp ; (sp+1):2*sp ; (2*sp+1):3*sp]; idx=idx';
        else
            %Layer potential f(r,n)/r
            SFr = [SV{1}(R,P) SV{2}(R,P) SW(R,P) SX(R,P)]./repmat(R+(R==0),1,4);
            %Pressure
            PrFr = FV{3}(R,P);
            MPr = PrFr; %.*(reshape(repmat(Vshc(1:sp,:),ntrg,1),sp,[]).');

            %indices
            idprin = 1:sp;
            idpr = sp+1:2*sp;
            idx = [1:sp ; (sp+1):2*sp ; (sp+1):2*sp ; (2*sp+1):3*sp]; idx=idx';
        end

        eu = repmat(sin(u),1,3).*[cos(u).*cos(v) cos(u).*sin(v) -sin(u)];
        ev = repmat(sin(u),1,3).*[-sin(v) cos(v) zeros(size(u,1),1)];

        Ndotr = sum(Nr(1:3:end,:).*er,2); Ndotr = repmat(Ndotr.',3,1);
        Ndotu = sum(Nr(1:3:end,:).*eu,2)./sin(u).^2;
        Ndotu = repmat(Ndotu.',3,1);
        Ndotv = sum(Nr(1:3:end,:).*ev,2)./sin(u).^2;
        Ndotv = repmat(Ndotv.',3,1);

        Ndote = [Ndotu(:) Ndotv(:) Ndotr(:)];

        % Radial and Tangential contributions
        MRqh = Fr; %.*(reshape(repmat(qhex,ntrg,1),4*sp,[]).');
        MGqh = SFr; %.*(reshape(repmat(qhex,ntrg,1),4*sp,[]).');

        Mqh = zeros(ntrg*d2,3*sp);
        Y = zeros(ntrg*d2,3*np);
        VWX = zeros(ntrg*d2,3*sp);
        for k=0:p
            ind = k^2+1:(k+1)^2;
            indv = [ind ind+sp ind+2*sp];
            Vk = Vnm('Vnm',Sc,k,[],u,v,er);
            Wk = Vnm('Wnm',Sc,k,[],u,v,er);
            Xk = Vnm('Xnm',Sc,k,[],u,v,er);
            VWX(:,indv) = [Vk Wk Xk];
        end

        for l=1:ntrg*d2
            for j=1:3
                Y(l,:) = Y(l,:) + Ndote(l,j)*( VWX(l,:)*( ...
                    ShTg{j}(:,idx(:))*( repmat(MGqh(l,:).',1,3*np).*T(idqh,:) ) + ...
                    ShTr{j}(:,idx(:))*( repmat(MRqh(l,:).',1,3*np).*T(idqh,:) ) + ...
                    ShY{j}(:,idpr)*(repmat(MPr(l,:).',1,3*np).*T(idprin,:)) ));
            end
        end

        Y = reshape(Y,ntrg,[]);
    end

    if size(Mqh,2)==3*sp
        %{
    Y = zeros(ntrg*d2,3*sp); 
    for k=0:p
        ind = k^2+1:(k+1)^2; 
        indv = [ind ind+sp ind+2*sp]; 
        Vk = Vnm('Vnm',Sc,k,[],u,v,er); 
        Wk = Vnm('Wnm',Sc,k,[],u,v,er); 
        Xk = Vnm('Xnm',Sc,k,[],u,v,er); 
   
        Yk = [Vk Wk Xk].*Mqh(:,indv);
        Y(:,indv) = Yk; 
    end
        %}
    else
        Y = zeros(ntrg*d2,4*sp);
        T = T(idqh,:);
        for k=0:p
            ind = k^2+1:(k+1)^2;
            indv = [ind ind+sp ind+2*sp ind+3*sp];
            Vk = Vnm('Vnm',Sc,k,[],u,v,er);
            Wk = Vnm('Wnm',Sc,k,[],u,v,er);
            Xk = Vnm('Xnm',Sc,k,[],u,v,er);

            if out
                Yk = [Vk Vk Wk Xk].*Mqh(:,indv);
            else
                Yk = [Vk Wk Wk Xk].*Mqh(:,indv);
            end

            Y(:,indv) = Yk;
        end

        tmp = Y; Y = zeros(ntrg,3*np*d2);
        for i=1:d2
            ind = (1:3*np) + (i-1)*3*np;
            if ~isempty(Rdiag)
                Y(:,ind) = reshape(Rdiag(ind(1))*tmp(ind,:)*T,ntrg,[]);
            else
                Y(:,ind) = tmp*T;
            end
        end

    end

end

end


%% Helper function to load kernel memory from file
function [ShTg,ShTr,ShY] = get_Surfgrad_coeffs(px,p)

printMsg('* Reading the Traction kernel coefficient matrices for p=%d: ',p);
[TrSz,rFlag] = readData('TractionCoeffsSpSize',p,1);

if(~rFlag)
    printMsg('Failed\n');
    printMsg('* Generating the Traction kernel coefficient matrices for p=%d\n',p);

    [ShTg,ShTr,ShY] = Surfgrad_coeffs(px,p,1,'VW','VW',0);

    TrData = zeros(0,5);
    for j=1:3
        [ig,jg,vg] = find(ShTg{j}); ng = [length(vg) ; zeros(length(vg)-1,1)];
        [ir,jr,vr] = find(ShTr{j}); nr = [length(vr) ; zeros(length(vr)-1,1)];
        [iy,jy,vy] = find(ShY{j}); ny = [length(vy) ; zeros(length(vy)-1,1)];
        TrData = [TrData ; [ig ; ir ; iy] [jg ; jr ; jy] [real(vg) ; real(vr) ; real(vy)] [imag(vg) ; imag(vr) ; imag(vy)] [ng ; nr ; ny]]; %#ok<AGROW>
    end

    writeData('TractionCoeffsSp',p, TrData);
    writeData('TractionCoeffsSpSize',p, size(TrData,1));
    printMsg('* Stored generated Traction coefficient matrices for p=%d\n', p);
else
    sp = (p+1)^2;
    [TrData,~] = readData('TractionCoeffsSp',p,[TrSz 5]);
    ShTg = cell(3,1); ShTr = ShTg; ShY = ShTg;
    lv = [find(TrData(:,5)>0) ; TrSz];
    sp3 = 3*sp;

    for j=1:3
        ind = lv(1+3*(j-1)):(lv(2+3*(j-1))-1);
        ShTg{j} = sparse(TrData(ind,1),TrData(ind,2),complex(TrData(ind,3),TrData(ind,4)),sp3,sp3);
        ind = lv(2+3*(j-1)):(lv(3+3*(j-1))-1);
        ShTr{j} = sparse(TrData(ind,1),TrData(ind,2),complex(TrData(ind,3),TrData(ind,4)),sp3,sp3);
        ind = lv(3+3*(j-1)):(lv(4+3*(j-1))-1);
        ShY{j} = sparse(TrData(ind,1),TrData(ind,2),complex(TrData(ind,3),TrData(ind,4)),sp3,sp3);
    end

    printMsg('Successful\n');
end


end

%% Helper funciton to obtain function handles of a specific type
function [FV,FW,FX] = LOCAL_get_Frn(type,out)
% Eigenvalues
SVeg  = @(n) n./((2*n+1).*(2*n+3));
SWeg  = @(n) (n+1)./((2*n+1).*(2*n-1));
SXeg  = @(n) 1./(2*n+1);
if out
    switch type
        case 'SMat'
            % f(r,n) exterior (SMat)
            SVext  = @(r,n) SVeg(n).*(r.^(-n-2));
            SWVext = @(r,n) (n./(4*n+2)).*(r.^(-n-2)-r.^(-n));
            SWWext = @(r,n) SWeg(n).*r.^(-n);
            SXext  = @(r,n) SXeg(n).*r.^(-n-1);

            FV = SVext;
            FW = {SWWext,SWVext};
            FX = SXext;
        case 'SpMat'
            % f'(r,n) exterior (SpMat)
            SpVext  = @(r,n) -(n+2).*SVeg(n).*(r.^(-n-3));
            SpWVext = @(r,n) (n./(4*n+2)).*(-(n+2).*r.^(-n-3)+n.*r.^(-n-1));
            SpWWext = @(r,n) -n.*SWeg(n).*r.^(-n-1);
            SpXext  = @(r,n) -(n+1).*SXeg(n).*r.^(-n-2);
            PWext = @(r,n) -n.*r.^(-n-1);

            FV = SpVext;
            FW = {SpWWext,SpWVext,PWext};
            FX = SpXext;

        case 'TMat'
            % f(r,n) exterior (Traction)
            %TV
            TVext = @(r,n) 2*(-n-2).*SVeg(n).*r.^(-n-3);
            %TW
            %Simplified these two as f(n)r^(-n-1)+g(n)r^(-n-3)
            TWVext = @(r,n) ((n.*(n+2))./(2*n+1)).*(-r.^(-n-3)+r.^(-n-1));
            TWWext = @(r,n) ((1+2*n.^2)./(1-4*n.^2)).*r.^(-n-1);

            %TX
            TXext = @(r,n) -(n+2).*SXeg(n).*r.^(-n-2);

            FV = TVext;
            FW = {TWWext,TWVext};
            FX = TXext;
        case 'DMat'
            % f(r,n) exterior (double layer)
            %DV
            DVext = @(r,n) ((2*n.^2+4*n+3)./((2*n+1).*(2*n+3))).*r.^(-n-2);
            %DW
            DWVext = @(r,n) ((n.*(n-1))./(2*n+1)).*(r.^(-n-2)-r.^(-n));
            DWWext = @(r,n) ((2*n.^2-2)./(4*n.^2-1)).*r.^(-n);

            %DX
            DXext = @(r,n) ((n-1)./(2*n+1)).*r.^(-n-1);

            FV = DVext;
            FW = {DWWext,DWVext};
            FX = DXext;
        case 'DpMat'
            % f(r,n) exterior (double layer)
            %DV
            DVext = @(r,n) (-n-2).*((2*n.^2+4*n+3)./((2*n+1).*(2*n+3))).*r.^(-n-3);
            %DW
            DWVext = @(r,n) ((n.*(n-1))./(2*n+1)).*((-n-2).*r.^(-n-3)+n.*r.^(-n-1));
            DWWext = @(r,n) -n.*((2*n.^2-2)./(4*n.^2-1)).*r.^(-n-1);

            %DX
            DXext = @(r,n) (-n-1).*((n-1)./(2*n+1)).*r.^(-n-2);

            %Pressure
            PWext = @(r,n) -n.*(n+1).*r.^(-n-2);

            FV = DVext;
            FW = {DWWext,DWVext,PWext};
            FX = DXext;
        case 'SDMat'
            % f(r,n) exterior (single + double layer)
            %SDV
            DVext = @(r,n) ((n+1)./(2*n+1)).*r.^(-n-2);
            %SDW
            DWVext = @(r,n) (n.*(2*n-1)./(4*n+2)).*(r.^(-n-2)-r.^(-n));
            DWWext = @(r,n) ((n+1)./(2*n+1)).*r.^(-n);

            %SDX
            DXext = @(r,n) (n./(2*n+1)).*r.^(-n-1);

            FV = DVext;
            FW = {DWWext,DWVext};
            FX = DXext;
    end
else

    switch type
        case 'SMat'
            % f(r,n) interior (SMat)
            SVVint = @(r,n) SVeg(n).*r.^(n+1);
            SVWint = @(r,n) -((n+1)./(4*n+2)).*(r.^(n-1)-r.^(n+1));
            SWint  = @(r,n) SWeg(n).*r.^(n-1);
            SXint  = @(r,n) SXeg(n).*r.^n;

            FV = {SVVint,SVWint};
            FW = SWint;
            FX = SXint;
        case 'SpMat'
            % f'(r,n) interior (SpMat)
            SpVVint = @(r,n) (n+1).*SVeg(n).*r.^n;
            SpVWint = @(r,n) -((n+1)./(4*n+2)).*((n-1).*r.^(n-2)-(n+1).*r.^n);
            SpWint  = @(r,n) (n-1).*SWeg(n).*r.^(n-2);
            SpXint  = @(r,n) n.*SXeg(n).*r.^(n-1);
            PVint = @(r,n) -(n+1).*r.^n;

            FV = {SpVVint,SpVWint,PVint};
            FW = SpWint;
            FX = SpXint;

        case 'TSMat'
            % f(r,n) interior (Traction of SL)
            %TV
            %Simplified these two as f(n)r^n+g(n)r^(n-2)
            TVVint = @(r,n) ((3+4*n+2*n.^2)./(3+8*n+4*n.^2)).*r.^n;
            TVWint = @(r,n) ((-1+n.^2)./(1 + 2*n)).*(r.^n-r.^(n-2));
            %TW
            TWint = @(r,n) 2*(n-1).*SWeg(n).*r.^(n-2);
            %TX
            TXint = @(r,n) (n-1).*SXeg(n).*r.^(n-1);

            FV = {TVVint,TVWint};
            FW = TWint;
            FX = TXint;
        case 'DMat'
            % f(r,n) interior (Double Layer)
            %DV
            DVVint = @(r,n) ((-2*n.*(n+2))./((2*n+1).*(2*n+3))).*r.^(n+1);
            DVWint = @(r,n) -(((n+1).*(n+2))./(2*n+1)).*(r.^(n+1)-r.^(n-1));
            %TW
            DWint = @(r,n) -((2*n.^2+1)./((2*n+1).*(2*n-1))).*r.^(n-1);
            %TX
            DXint = @(r,n) -((n+2)./(2*n+1)).*r.^n;

            FV = {DVVint,DVWint};
            FW = DWint;
            FX = DXint;
        case 'DpMat'
            % f(r,n) interior (normal derivative Double Layer)
            %DV
            DVVint = @(r,n) (n+1).*((-2*n.*(n+2))./((2*n+1).*(2*n+3))).*r.^n;
            DVWint = @(r,n) -(((n+1).*(n+2))./(2*n+1)).*((n+1).*r.^n-(n-1).*r.^(n-2));
            %TW
            DWint = @(r,n) -(n-1).*((2*n.^2+1)./((2*n+1).*(2*n-1))).*r.^(n-2);
            %TX
            DXint = @(r,n) -n.*((n+2)./(2*n+1)).*r.^(n-1);

            PVint = @(r,n) n.*(n+1)*r.^(n-1);

            FV = {DVVint,DVWint,PVint};
            FW = DWint;
            FX = DXint;
        case 'SDMat'
            % f(r,n) interior (SL + DL)
            %[S+D]V
            DVVint = @(r,n) (-n./(2*n+1)).*r.^(n+1);
            DVWint = @(r,n) (-(n+1).*(2*n+3)./(4*n+2)).*(r.^(n+1)-r.^(n-1));
            %[S+D]W
            DWint = @(r,n) (-n./(2*n+1)).*r.^(n-1);
            %[S+D]X
            DXint = @(r,n) (-(n+1)./(2*n+1)).*r.^n;

            FV = {DVVint,DVWint};
            FW = DWint;
            FX = DXint;
    end
end

end
