pstart=2; pend=16; 
parr=pstart:pend; spend=(pend+1)^2;
uu = 2/sqrt(3);
a=1./uu;

%------------------------------------------------------------------------%

DY_cell=cell(1,length(parr));
SY_cell=cell(1,length(parr));
DYsh_cell=cell(1,length(parr));
SYsh_cell=cell(1,length(parr));
myDY_cell=cell(1,length(parr));
mySY_cell=cell(1,length(parr));
myDYsh_cell=cell(1,length(parr));
mySYsh_cell=cell(1,length(parr));

for l=1:length(parr)
    
    p=parr(l);

    %----------------------------------------------------------%
    % Compute spheroidal harmonic functions and eigenvalues from formula
    %----------------------------------------------------------%
    np = 2*p*(p+1); sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    [u,v]=gl_grid(p);
    
    Y = zeros(np,sp); SYns = Y;  Y2=Y; 
    myDY=zeros(np,sp); mySY=zeros(np,sp);
    DLeigval=zeros(sp,1);
    SLeigval=zeros(sp,1);

    %Compute spectra using formulas
    [~,DLeigval,~]=DLspectrum(p,u0);
    [~,SLeigval,~]=SLspectrum(p,u0,a);
    
    for k=1:sp
        n = nn(k); m = mm(k); 
        Y(:,k) = Ynm(n,m,u,v); 
        fac=(1./sqrt(uu^2-cos(u).^2)); fac=fac(:); 
        Y2(:,k) = fac.*Y(:,k);
    end

    %Form exact S[Y2] and D[Y] matrices using formulas
    myDY=repmat(DLeigval',np,1).*Y;
    mySY=repmat(SLeigval',np,1).*Y2;
    
    %----------------------------------------------------------%
    % Compute S[Y2] and D[Y], along with spheroidal harmonic 
    % coefficients, using singular quadrature

    [DYsh,DY]=SDY(p,'DL');
    [SYsh,SY]=SDY(p,'SL');

    %----------------------------------------------------------%
    %store eigenvalues in a diagonal matrix
    myDYsh=diag(DLeigval);
    mySYsh=diag(SLeigval);
    
    %----------------------------------------------------------%
    %Store all matrices in cell

    DY_cell{l}=DY;
    SY_cell{l}=SY;
    DYsh_cell{l}=DYsh;
    SYsh_cell{l}=SYsh;
    myDY_cell{l}=myDY;
    mySY_cell{l}=mySY;
    myDYsh_cell{l}=myDYsh;
    mySYsh_cell{l}=mySYsh;

end

%%

% PLOTTING
%------------------------------------------------------------%

close all 

DYerr_cell=cell(1,length(parr));
SYerr_cell=cell(1,length(parr));
DYsherr_cell=cell(1,length(parr));
SYsherr_cell=cell(1,length(parr));

% Seeing how the error changes with order p, for all eigs
[d,~]=size(DYsherr_cell{end});
Dconv=zeros(d,length(parr));
Sconv=zeros(d,length(parr));

Sratio_cell=cell(1,length(parr));

for l=1:length(parr)
    DYerr_cell{l}=abs(DY_cell{l}-myDY_cell{l});
    SYerr_cell{l}=abs(SY_cell{l}-mySY_cell{l});
    DYsherr_cell{l}=abs(DYsh_cell{l}-myDYsh_cell{l});
    SYsherr_cell{l}=abs(SYsh_cell{l}-mySYsh_cell{l});

    [sl,~]=size(SYsherr_cell{l});
    mySYsh_l=mySYsh_cell{l};
    SYsh_l=SYsh_cell{l};
    Sratio_l=zeros(sl);
    for i=1:sl
        for j=1:sl
            if SYsh_l(i,j) ~= 0
                Sratio_l(i,j)=mySYsh_l(i,j)/SYsh_l(i,j);
            end
        end
    end
    Sratio_cell{l}=Sratio_l;

    % Seeing how the error changes with order p, for all eigs
    Deigs_l=diag(DYsherr_cell{l});
    Dconv(1:length(Deigs_l),l)=Deigs_l;
    Seigs_l=diag(SYsherr_cell{l});
    Sconv(1:length(Seigs_l),l)=Seigs_l;


    if parr(l)==pend

    figure;
    imagesc(log10(DYerr_cell{l}))
    title('D[Y] error')
    colorbar
    ylabel("(n,m) index")
    xlabel("order of spheroidal harmonic expansion")
    
    figure;
    imagesc(log10(DYsherr_cell{l}))
    title('D[Y] coefficients error')
    colorbar
    ylabel("(n,m) index")
    xlabel("order of spheroidal harmonic expansion")

    figure;
    imagesc(log10(SYerr_cell{l}))
    title('S[Y2] error')
    colorbar
    ylabel("(n,m) index")
    xlabel("order of spheroidal harmonic expansion")

    figure;
    imagesc(log10(SYsherr_cell{l}))
    title('S[Y2] coefficients error')
    colorbar
    ylabel("(n,m) index")
    xlabel("order of spheroidal harmonic expansion")

    end
end

figure;
imagesc(log10(Dconv),x=parr)
colorbar
ylabel("(n,m) index")
xlabel("order of spheroidal harmonic expansion")
title("double layer eigenvalue error")

figure;
imagesc(log10(Sconv),x=parr)
colorbar
ylabel("(n,m) index")
xlabel("order of spheroidal harmonic expansion")
title("single layer eigenvalue error")

% figure;
% imagesc(abs(Sratio_l)); 
% colorbar;
% title("Ynm coefficients from formula / Ynm coefficients from singular quadrature")
% 
% figure;
% plot(1:sl,diag(real(Sratio_l)));
% title("Eigenvalues from formula / Eigenvalues from singular quadrature")
% xlabel("(n,m) index of matrix")
% ylabel("Re(ratio of Single layer eigenvalues)")
% 
% figure;
% plot(1:sl,diag(imag(Sratio_l)));
% title("Imaginary")
% 
% 
% Sratio_12=Sratio_cell{11};
% [d12,~]=size(Sratio_12);
% figure;
% plot(1:d12,diag(real(Sratio_12)));
% title("Eigenvalues from formula / Eigenvalues from singular quadrature, p=12")
% xlabel("(n,m) index of matrix")
% ylabel("Re(ratio of Single layer eigenvalues)")


