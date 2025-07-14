function PQcell = legendre_otc(p,u,Qoption,dPoption,dQoption,debug)
%computes the associated Legendre functions up to degree p and order 
% m =-p:p, evaluated for u>1.
%     (o) Qoption: bool that tells whether to compute Q_n^m
%     (o) dPoption: bool that tells whether to compute first derivs P_n^m
%     (o) dQoption: bool that tells whether to compute first derivs Q_n^m.
%         Must compute Q_n^m if dQoption=true.
%     (o) debug: prints intermediate results for debugging.

if nargin==5
    %optional input, if debug=1 it will print intermediate results.
    debug=0;
elseif nargin==4
    %optional input, can be set to zero if only Pnm is needed. By default
    %we also compute Qnm
    dQoption=0;
    debug=0;
elseif nargin==3
    dPoption=0;
    dQoption=0;
    debug=0;
elseif nargin==2
    Qoption=1;
    dPoption=0;
    dQoption=0;
    debug=0;
end
if p<0
    error('invalid degree; p must be nonnegative')
end
if Qoption==0 && dQoption==1
    Qoption=1;
    fprintf('Q option = 0 overriden to compute Q derivatives...\n')
end
u=reshape(u,1,length(u));
if all(imag(u)==0)
    if sum(u<=1)>0
        error('invalid u; u must be >1')
    end
end


PQcell=cell(1,Qoption+dPoption+dQoption+1); %cell to store P, Q, dP, dQ matrices

if dPoption || dQoption
    pmax=p+1;   %To compute derivatives we need to go an extra order
else
    pmax=p;
end
sp=(pmax+1)^2; %number of legendre fxns for order-p expansion
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
geti= @(n,m) m+n^2+n+1; % Map (n,m) to 0 <= k <= sp


% Calculate P_0^0 from simple formulas
%-----------------------------------------------------------------%
P=zeros(sp,length(u));
P(1,:)=ones(1,length(u)); %P00

if debug
fprintf("Computed P_0^0 \n")
if length(u)<5
fprintf("P_0^0(u)=\n")
disp(P(1,:))
end
end

%populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
%-----------------------------------------------------------------%
if pmax>0
    if debug
    fprintf('Calculating P_m^m for m=0:%d and P_{m+1}^m for m=0:%d... \n',p,p-1)
    fprintf('---------------------------\n')
    end

    for m=0:pmax

        if debug
        fprintf("this m=%d\n",m)
        end
        
        %P_m^m(u)=(2m-1)!! (u^2-1)^(m/2) if m>0
        if m>0
            Pmm=prod(1:2:2*m-1).*(u.^2-1).^(m/2); %P_m^m
            P(geti(m,m),:)=Pmm; 
            P(geti(m,-m),:)=1./factorial(2*m).*Pmm; %P_m^{-m}
            
            %testing/debug
            if debug
            fprintf("Computed P_%d^%d\n",[m,m])
            if length(u)<5
            fprintf("P_%d^%d(u)= \n",[m,m])
            disp(Pmm)
            end
            end  
        end

        %Formula: P_{m+1}^m = (2m+1)u P_m^m(u) 
        Pmp1m = (2*m+1).* u.* P(geti(m,m),:);

        if m<pmax
            P(geti(m+1,m),:)=Pmp1m;
            P(geti(m+1,-m),:)=1./factorial(2*m+1).*Pmp1m;

            %testing/debug
            if debug
            fprintf("computed P_%d^%d\n",[m+1,m])
            if length(u)<5       
            fprintf("P_%d^%d(u)= \n",[m+1,m])
            disp(Pmp1m)
            end
            end
        end
    end
    if debug; fprintf('---------------------------\n'); end
end
%-----------------------------------------------------------------%

%Use recursion if higher order is needed.
if pmax>1   
    
    if debug
    fprintf('Computing P_n^m for m=0:%d, n=m+2:%d \n',pmax-2,pmax)
    fprintf('---------------------------\n')
    end

    for m=0:pmax-2
    
        if debug
        fprintf('this m=%d\n',m)
        end

        %Use forward recursion to compute more P_n^m.
        %----------------------------------------------------------------%
        for n=m+2:pmax
       
            Pnm=( (2*n-1).*u.*P(geti(n-1,m),:)-(n+m-1).*P(geti(n-2,m),:))./(n-m);
            P(geti(n,m),:)=Pnm;
            P(geti(n,-m),:)=factorial(n-m)./factorial(n+m).*Pnm;

            if debug
            fprintf("Computed P_%d^%d\n",[n,m])
            if length(u)<5
            fprintf("P_%d^%d(u)= \n",[n,m])
            disp(Pnm)
            end
            end
        end    
        %-----------------------------------------------------------------
    end
    if debug
    fprintf('---------------------------\n')
    end
end
    
% Store P matrix into output cell 
PQcell{1}=P(1:(p+1)^2,:);

%If desired, calculate Q
%---------------------------------------------------------------------%
if Qoption
    Q=zeros(size(P));

    % Calculate highest order separately
    Ppp1p=(2*pmax+1).* u.*P(geti(pmax,pmax),:); %This will be needed to get Q_p^p
    Hp=cf(pmax+1,pmax,u);
    Qpp=factorial(2*pmax).*(-1).^pmax./(Ppp1p-Hp.*P(geti(pmax,pmax),:));
    Q(geti(pmax,pmax),:)=Qpp;
    Q(geti(pmax,-pmax),:)=1./factorial(2*pmax).*Qpp;

    if pmax>0
        for m=0:pmax-1
            Ppm=P(geti(pmax,m),:);
            Ppm1m=P(geti(pmax-1,m),:);
    
            H=cf(pmax,m,u);
            Qpm1m=factorial(pmax+m-1)/factorial(pmax-m).*(-1)^m./(Ppm-H.*Ppm1m);
            Qpm=H.*Qpm1m;
    
            Q(geti(pmax-1,m),:)=Qpm1m;
            Q(geti(pmax,m),:)=Qpm;

            Q(geti(pmax-1,-m),:)=factorial(pmax-1-m)./factorial(pmax-1+m).*Qpm1m;
            Q(geti(pmax,-m),:)=factorial(pmax-m)./factorial(pmax+m).*Qpm;
            
            if pmax>1
                % Use backward recursion to compute Q_n^m for n=p-2:m
                for n=pmax-2:-1:m
                    Qnm=((2*n+3).*u.*Q(geti(n+1,m),:)-(n-m+2).*Q(geti(n+2,m),:))./(n+m+1);
                    Q(geti(n,m),:)=Qnm;
                    Q(geti(n,-m),:)=factorial(n-m)./factorial(n+m).*Qnm;
            
                    if debug
                    fprintf("Computed Q_%d^%d \n",[n,m])
                    if length(u)<5
                    fprintf("Q_%d^%d(u)=\n",[n,m])
                    disp(Qnm)
                    end
                    end
                end
            end    
        end
    end

    PQcell{2}=Q(1:(p+1)^2,:);
end
%---------------------------------------------------------------------%

% Compute first derivatives
%---------------------------------------------------------------------%
if dPoption | dQoption
    if dPoption
        dP=zeros((p+1)^2,length(u));
    end
    if dQoption
        dQ=zeros((p+1)^2,length(u));
    end
    for k=1:(p+1)^2
        n=nn(k); m=mm(k);
    
        if dPoption
            dPnm=((m-n-1).*P(geti(n+1,m),:)+(n+1).*u.*P(geti(n,m),:))./(1-u.^2);
            dP(geti(n,m),:)=dPnm;

            if debug
            fprintf('Computed dP_%d^%d\n',n,m)
            end
        end

        if dQoption
            dQnm=((m-n-1).*Q(geti(n+1,m),:)+(n+1).*u.*Q(geti(n,m),:))./(1-u.^2);
            dQ(geti(n,m),:)=dQnm;

            if debug
            fprintf('Computed dQ_%d^%d\n',n,m)
            end
        end
    end
    if dPoption
        PQcell{3}=dP;
    end
    if dQoption
        PQcell{4}=dQ; %if calculated dP in index 3, put dQ in index 4
    end
end
%---------------------------------------------------------------------%


% Compute second derivatives.
% Basically the same code as above, but below is written with the purpose of
% backwards compatibility with the rest of the codebase.
% Note that we shall need access to first derivatives for the calculation
% of the second derivatives.
if dPoption == 2 | dQoption == 2
    if dPoption == 2
        d2P=zeros((p+1)^2,length(u));
    end
    if dQoption == 2
        d2Q=zeros((p+1)^2,length(u));
    end

    for k=1:(p+1)^2
        n = nn(k); m = mm(k);

        if dPoption == 2
            one_minus_x2 = 1 - u.^2;
            ddPnm_term = 2.*u.*dP(geti(n,m),:) - (n*(n + 1) - m^2./one_minus_x2).*P(geti(n,m),:);
            ddPnm = ddPnm_term ./ one_minus_x2;
            d2P(geti(n,m),:) = ddPnm;
        end

        if dPoption == 2
            % Alternate evaluation
        end

        if dQoption == 2
            one_minus_x2 = 1 - u.^2;
            ddQnm_term = 2.*u.*dQ(geti(n,m),:) - (n*(n + 1) - m^2./one_minus_x2).*Q(geti(n,m),:);
            ddQnm = ddQnm_term ./ one_minus_x2;
            d2Q(geti(n,m),:) = ddQnm;
        end
    end

    % Add them to the resulting cell (indices 5 and 6)
    if dPoption == 2
        PQcell{5} = d2P;
    end

    if dQoption == 2
        PQcell{6} = d2Q;
    end
end

%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
% Continued fraction functions

function [H,iters]=cf(n,m,x)
    %Compute the continued fraction that relates Q_n^m with Q_{n-1}^m, to
    %desired tolerance, using Modified Lentz's Method.

    tol=1e-15;
    % max_iters=100000;
    max_iters=1e5;


    tiny=1e-300;
    f0=tiny;
    c0=f0;
    d0=0;

    %sprintf('Initial values:')
    %sprintf("f_0=%.6e",f0)
    %sprintf("c_0=%.6e",c0)
    %sprintf("d_0=%.6e",d0)

    for k=1:max_iters
        %sprintf("iteration %d",k-n+1)

        d1=b(n,k,m,x)+a(n,k,m).*d0;
        if d1==0
            d1=tiny;
        end

        c1=b(n,k,m,x)+a(n,k,m)./c0;
        if c1==0
            c1=tiny;
        end

    %   sprintf("c_%d=%.6e",k-n+1,c1)
    %   sprintf("d_%d=%.6e",k-n+1,d1)

        d1=1./d1;
        delta=c1.*d1;
        f1=delta.*f0;

    %   sprintf("f_%d=%.6e",k-n+1,f1)
    %   sprintf("abs(1-delta)=%.6e",abs(1-delta))

        iters=k;

        if abs(delta-1)<tol
            break
        end

        if k==max_iters
            error("Continued fraction algorithm did not converge.")
        end

        %update values for recursion
        c0=c1;
        d0=d1;
        f0=f1;
    
    end
    H=f1;
end

function bj=b(n,j,m,x)
    %b_j(n,m,x) coefficient in recursion relation
    % bj= (2*(n+j-1)+1).*x ./(n+j-1+m);

    bj=-(2*(n+j-1)+1).*x ./(n+j-m);
end

function aj=a(n,j,m)
    %a_j coefficient in recursion relation
    % if j==1
    %     aj=1;
    % else
    %     aj=-(n+j-1-m)./(n+j-2+m);
    % end

    aj=-(n+j-1+m)./(n+j-m);
end

end