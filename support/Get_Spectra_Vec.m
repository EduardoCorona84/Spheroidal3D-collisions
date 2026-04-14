function [Spectra] = Get_Spectra_Vec(p,lambda,r,type,out)
%               Computed Spherical Harmonic coefficients
%inputs:        p (int>0): degree of spherical harmonics
%               lambda (double>0): parameter in modified Laplace
%               r (double nr x 1):  distance to target points
%               type (string): SMat/DMat/SpMat/DpMat-which kernel
%               out (bool): 1 if exterior, 0 if interior 
%outputs:       Spectra (double nr x p+1): SpHarm coefficients for degree
%                   0-p at each entry in r

r=reshape(r,1,[]);
%Calls subfunction depending on type
if strcmp(type,'SMat')
    Spectra=Spectra_nr2_local(p,lambda,r,out);
elseif strcmp(type,'DMat')
    if out
        Spectra=SpectraDLO_nr2_local(p,lambda,r);
    else
        Spectra=SpectraDLI_nr_local(p,lambda,r);
    end
elseif strcmp(type,'SDMat')
    Spectra=Spectra_nr2_local(p,lambda,r,out);
    display(Spectra)
    if out
        Spectra=Spectra+SpectraDLO_nr2_local(p,lambda,r);
    else
        Spectra=Spectra+SpectraDLI_nr_local(p,lambda,r);
        display(Spectra)
    end
elseif strcmp(type,'SpMat')
        Spectra=SpectraSP_local(p,lambda,r,out);
elseif strcmp(type,'DpMat')
        Spectra=SpectraDP_local(p,lambda,r,out);
end
end


function [SL]=Spectra_nr2_local(n,lambda,r,out)
%Computes Single Layer Spectra
%inputs: n:p from main
%        lambda: from main
%        r: from main
%        out: from main
% Output: single layer spectra 
lr=lambda*r;
p = (0:n)+0.5;
dm=n;
ilambda = (1i)*lambda;
SL=zeros(size(r,2),n+1);
mxlr=5; 
ind=lr<mxlr;
rless=r(ind);
rmore=r(~ind);

if nargin<4
    out=1; 
end


if out
    
    % modified Bessel functions of first and second kind
    %J_n(lambda)
    J_pos = besselj(p,ilambda);
    MJ_first_pos = (1i.^(-p)).*J_pos;
    
    
    MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
%*************************************************************************
%*************************************************************************
    %uses Matlab functions for 2nd sph bessel if r*lambda< maximum  
    if ~isempty(rless) 
        nr=size(rless,2);
        Rless=repmat(rless',1,n+1);
        P=repmat(p,nr,1);
     
        Jr_pos=besselj(P,ilambda*Rless);
        Jr_neg=besselj(-P,ilambda*Rless);
        
        MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
        MJr_first_neg =bsxfun(@times,Jr_neg,(1i.^p));
        MJr_second = bsxfun(@times,(MJr_first_neg - MJr_first_pos),(pi./(2*sin(p*pi))));
        
        % modified spherical Bessel functions of first and second kind
        MSJ_first = MJ_first_pos; %sqrt(pi/(2*lambda))*MJ_first_pos;
        MSJr_second = bsxfun(@times,sqrt(1./rless)',MJr_second); %sqrt(pi/(2*lambda*r))*MJr_second;
        
        %spectra
        SL(ind,:) = real(bsxfun(@times,MSJr_second,MSJ_first));
    end
%uses asymptotic approximation if %r*lambda larger than maximum     
%*************************************************************************
%*************************************************************************    
    if ~isempty(rmore)
       nr=size(rmore,2);
        %approximation k_n(z) ~ (pi/2z)*exp(-z) as z->infty
        MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
        MSJr_second = bsxfun(@times,(exp(-lambda*rmore)./rmore)',ones(nr,n+1));
       
        LRMore=repmat((1./(lambda*rmore))',1,n+1);
        EXP=repmat((0:dm),nr,1);
        LR=LRMore.^EXP;
        
        
        ank = zeros(n+1,dm+1); 
        ank(1,:) = [1 zeros(1,dm)]; 
        
        for j=1:n
            nck = zeros(dm,1); 
            for i=0:min(j,dm)
                nck(i+1) = nchoosek(j,i);
                ank(j+1,i+1)=nck(i+1)*prod(j+1:j+i)/(2^i);
            end   
        end
        
        for i=1:nr
            SLRank(:,i) = sum(repmat(LR(i,:),n+1,1).*ank,2)';
        end
        
        MSJr_second = bsxfun(@times,SLRank,(exp(-lambda*rmore)./rmore))'; 
        
        SL(~ind,:) = real(bsxfun(@times,MSJr_second,MSJ_first)); 
    end    
   %}
    SL=real(bsxfun(@times,MSJr_second,MSJ_first));

%*************************************************************************
%*************************************************************************
else %The Interior Case 

    
    
    nr=size(r,2);
    R=repmat(r',1,n+1);
    P=repmat(p,nr,1);
    Jr_pos=besselj(P,ilambda*R);
    
    
    MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
    
    %K_n(lambda)
    J_pos = besselj(p,ilambda);
    J_neg = besselj(-p,ilambda);
    
    MJ_first_pos = (1i.^(-p)).*J_pos;
    MJ_first_neg = (1i.^p).*J_neg;
    
    MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);
    MJ_second=repmat(MJ_second,nr,1);
   
    % modified spherical Bessel functions of first and second kind
    MSJr_first = MJr_first_pos; %sqrt(pi/(2*lambda))*MJr_first_pos;
    MSJ_second = bsxfun(@times,sqrt(1./r)',MJ_second); %sqrt(pi/(2*lambda*r))*MJ_second;
    %spectra
    SL = real(MSJr_first.*MSJ_second);
    
end
end

function [DL_I,der_MSJ_second] = SpectraDLI_nr_local(n,lambda,r )
%Computes Double Layer Spectra in the Interior
%inputs: n:p from main
%        lambda: from main
%        r: from main
%        out: from main
% Output: Double layer spectra Interior
p = [0:n]+0.5;
ilambda = (1i)*lambda;

% modified Bessel functions of first and second kind
nr=size(r,2);
R=repmat(r',1,n+1);
P=repmat(p,nr,1);

Jr_pos=besselj(P,ilambda*R);
Jr_neg=besselj(-P,ilambda*R)   

nr=size(r,2);
R=repmat(r',1,n+1);
P=repmat(p,nr,1);

Jr_pos=besselj(P,ilambda*R);
MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));

J_pos = besselj(p,ilambda);
J_neg = besselj(-p,ilambda);
MJ_first_pos = (1i.^(-p)).*J_pos;
MJ_first_neg = (1i.^p).*J_neg;
MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);


MSJr_first = bsxfun(@times,sqrt(pi./(2*lambda*r))',MJr_first_pos);
MSJ_second = sqrt(pi/(2*lambda))*MJ_second;
der_MSJ_second = der_MSJ_s(MSJ_second, lambda);

% spectra of double layer (inner) 
DL_I = -real(bsxfun(@times,((-2*lambda^2*der_MSJ_second)/pi),MSJr_first));


end

function [ DLO ] =SpectraDLO_nr2_local( n,lambda,r )
%Computes Double Layer Spectra in the Exterior
%inputs: n:p from main
%        lambda: from main
%        r: from main
%        out: from main
% Output: Double layer spectra Exterior 
mxlr=5;
lr=lambda*r;
dm=n;
p = (0:n)+0.5;
ilambda = (1i)*lambda;
rlambda=lambda*r;
ind=rlambda<mxlr;

rless=r(ind);
rmore=r(~ind);

DLO=zeros(size(r,2),n+1);
J_pos = besselj(p,ilambda);

nr=size(r,2);
R=repmat(r',1,n+1);
P=repmat(p,nr,1);

Jr_pos=besselj(P,ilambda*R);
Jr_pos=besselj(P,ilambda*R);
Jr_neg=besselj(-P,ilambda*R);
J_neg = besselj(-p,ilambda);

MJ_first_pos = (1i.^(-p)).*J_pos;
MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
MJ_first_neg = (1i.^p).*J_neg;
MJr_first_neg =bsxfun(@times,Jr_neg,(1i.^p)); 

MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
der_MSJ_first = der_MSJ_f(MSJ_first, lambda);

%*************************************************************************
%*************************************************************************
%Use MatLab's function if lambda*r < maximum allowed
if ~isempty(rless)

    MJr_second = bsxfun(@times,(MJr_first_neg(ind,:) - MJr_first_pos(ind,:)),(pi./(2*sin(p*pi))));
    MSJr_second = bsxfun(@times,sqrt(pi./(2*lambda*rless))',MJr_second);
    DLO(ind,:) = -real(bsxfun(@times,MSJr_second,(-(2*lambda^2/pi)*der_MSJ_first)));

end
%*************************************************************************
%*************************************************************************
%Use asympotic approximation otherwise
if ~isempty(rmore)
    nr=size(rmore,2);
    LRMORE=repmat((1./(lambda*rmore))',1,n+1);
    EXP=repmat((0:dm),nr,1);
    LR=LRMORE.^EXP;
      
    ank = zeros(n+1,dm+1); 
    ank(1,:) = [1 zeros(1,dm)]; 
        
    for j=1:n
        nck = zeros(dm,1); 
        for i=0:min(j,dm)
            nck(i+1) = nchoosek(j,i);
            ank(j+1,i+1)=nck(i+1)*prod(j+1:j+i)/(2^i);
        end   
    end
    
    for i=1:nr
        SLRank(:,i) = sum(repmat(LR(i,:),n+1,1).*ank,2)';
    end
     
    MSJr_second = bsxfun(@times,SLRank,(exp(-lambda*rmore)./rmore))';       
    DLO(~ind,:) = -real(bsxfun(@times,MSJr_second,(-(lambda)*der_MSJ_first)));
end
end

function [SP]=SpectraSP_local(n,lambda,r,out)
%Computes S' Spectra 
%inputs: n:p from main
%        lambda: from main
%        r: from main
%        out: from main
% Output: S' Spectra

if out==1
    p = [0:n]+0.5;
    ilambda = (1i)*lambda;
    nr=size(r,2);
    R=repmat(r',1,n+1);
    P=repmat(p,nr,1);
    
    J_pos = besselj(p,ilambda);
    J_neg = besselj(-p,ilambda);
    
    MJ_first_pos = (1i.^(-p)).*J_pos;
    Jr_pos=besselj(P,ilambda*R);
    Jr_neg=besselj(-P,ilambda*R);
    MJ_first_neg = (1i.^p).*J_neg;
    MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
    MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
    MJr_first_neg = bsxfun(@times,Jr_neg,(1i.^p));
    
    MJR_second=repmat((pi./(2*sin(p*pi))),size(r,2),1).*(MJr_first_neg-MJr_first_pos);
    MSJR_second=sqrt(pi./(2*lambda*R)).*MJR_second;
    MSJR_second=local_second_correction(r,n,lambda);
    
    der_MSJR_second=der_MSJ_s(MSJR_second,lambda*r);
    SP=der_MSJR_second.*repmat(MSJ_first,size(der_MSJR_second,1),1)*(2*lambda^2)/pi;
elseif out==0
    p = [0:n]+0.5;
    ilambda = (1i)*lambda;
    nr=size(r,2);
    R=repmat(r',1,n+1);
    P=repmat(p,nr,1);
    
    Jr_pos=besselj(P,ilambda*R);
    Jr_neg=besselj(-P,ilambda*R);
    J_pos = besselj(p,ilambda);
    J_neg = besselj(-p,ilambda);
    
    MJ_first_pos = (1i.^(-p)).*J_pos;
    MJ_first_neg = (1i.^p).*J_neg;
    MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);
    MSJ_second = sqrt(pi/(2*lambda))*MJ_second;
    
    MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
    MSJr_first = bsxfun(@times,sqrt(pi./(2*lambda*r))',MJr_first_pos);
    der_MSJR_first=der_MSJ_f(MSJr_first,lambda*r);
    SP=der_MSJR_first.*repmat(MSJ_second,size(der_MSJR_first,1),1)*(2*lambda^2)/pi;    
end
end
function [DP]=SpectraDP_local(n,lambda,r,out)
%Computes D' Spectra 
%inputs: n:p from main
%        lambda: from main
%        r: from main
%        out: from main
% Output: D' Spectra
p = (0:n)+0.5;
ilambda = (1i)*lambda;
if out == 1
    der_MSJR_second=der_MSJ_s(local_second_correction(r,n,lambda),lambda*r);
    J_pos = besselj(p,ilambda);
    MJ_first_pos = (1i.^(-p)).*J_pos;
    MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
    der_MSJ_first=der_MSJ_f(MSJ_first,lambda);
    DP=(2*lambda^3/pi).*der_MSJR_second.*repmat(der_MSJ_first,size(r,2),1);
elseif out ==0
    nr=size(r,2);
    J_pos = besselj(p,ilambda);
    J_neg = besselj(-p,ilambda);
    MJ_first_pos = (1i.^(-p)).*J_pos;
    MJ_first_neg = (1i.^p).*J_neg;
    MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);
    MSJ_second = sqrt(pi/(2*lambda))*MJ_second;
    der_MSJ_second=der_MSJ_s(MSJ_second,lambda);
    R=repmat(r',1,n+1);
    P=repmat(p,nr,1);
    Jr_pos=besselj(P,ilambda*R);
    MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
    MSJr_first = bsxfun(@times,sqrt(pi./(2*lambda*r))',MJr_first_pos);
    der_MSJR_first=der_MSJ_f(MSJr_first,lambda*r);
    DP=(2*lambda^3/pi)*der_MSJR_first.*repmat(der_MSJ_second,nr,1);
end
end

function [MSJr_second]=local_second_correction(r,n,lambda)
%This funciton makes the correction for modified 2nd bessel function used
%in the D' case
%inputs: From SpectraDP_local
%Outputs: Modified 2nd Spherical Bessel Function

MSJR_second=zeros(size(r,2),(n+1));
mxlr=5;
lr=lambda*r;
dm=n;
if nargin<4
    out=1; 
end

p = (0:n)+0.5;
ilambda = (1i)*lambda;
rlambda=lambda*r;
ind=rlambda<mxlr;
rless=r(ind);
rmore=r(~ind);
J_pos = besselj(p,ilambda);
nr=size(r,2);
R=repmat(r',1,n+1);
P=repmat(p,nr,1);

Jr_pos=besselj(P,ilambda*R);
Jr_pos=besselj(P,ilambda*R);
Jr_neg=besselj(-P,ilambda*R);
J_neg = besselj(-p,ilambda);

MJ_first_pos = (1i.^(-p)).*J_pos;
MJr_first_pos = bsxfun(@times,Jr_pos,(1i.^(-p)));
MJ_first_neg = (1i.^p).*J_neg;
MJr_first_neg =bsxfun(@times,Jr_neg,(1i.^p)); 
MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;

der_MSJ_first = der_MSJ_f(MSJ_first, lambda);

%*************************************************************************
%*************************************************************************
%use Matlab's function if less than max lambda*r...
if ~isempty(rless)

    MJr_second = bsxfun(@times,(MJr_first_neg(ind,:) - MJr_first_pos(ind,:)),(pi./(2*sin(p*pi))));
    MSJr_second(ind,:) = bsxfun(@times,sqrt(pi./(2*lambda*rless))',MJr_second); 
end

%*************************************************************************
%*************************************************************************
%Otherwise, use asymptotic approximation
if ~isempty(rmore)
    nr=size(rmore,2);
    LRMORE=repmat((1./(lambda*rmore))',1,n+1);
    EXP=repmat((0:dm),nr,1);
    LR=LRMORE.^EXP;
    ank = zeros(n+1,dm+1); 
    ank(1,:) = [1 zeros(1,dm)]; 
        
    for j=1:n
        nck = zeros(dm,1); 
        for i=0:min(j,dm)
            nck(i+1) = nchoosek(j,i);
            ank(j+1,i+1)=nck(i+1)*prod(j+1:j+i)/(2^i);
        end   
    end
    for i=1:nr
        SLRank(:,i) = sum(repmat(LR(i,:),n+1,1).*ank,2)';
    end
     MSJr_second(~ind,:) = bsxfun(@times,SLRank,(exp(-lambda*rmore)./rmore))';       
end
end




