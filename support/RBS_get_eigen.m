function [FV,FW,FX] = RBS_get_eigen(type,out)
% Eigenvalues
%SVeg  = @(n) n./((2*n+1).*(2*n+3)); 
%SWeg  = @(n) (n+1)./((2*n+1).*(2*n-1));
%SXeg  = @(n) 1./(2*n+1);
if out
switch type 
    case 'SMat'
    % f(r,n) exterior (SMat) 
    FV = @(n) n./((2*n+1).*(2*n+3)); %SVeg; 
    FW = @(n) (n+1)./((2*n+1).*(2*n-1)); %SWeg; 
    FX = @(n) 1./(2*n+1); %SXeg; 
    case 'SpMat'
    % f'(r,n) exterior (SpMat)
    SpVext  = @(n) -(n+2).*(n./((2*n+1).*(2*n+3))); %-(n+2)*SVeg(n); 
    SpWVext = @(n) (-2*n./(4*n+2));
    SpWWext = @(n) -n.*((n+1)./((2*n+1).*(2*n-1))); %-n*SWeg(n);
    SpXext  = @(n) -(n+1)./(2*n+1); %-(n+1)*SXeg(n);
    PWext = @(n) -n; 

    FV = SpVext; 
    FW = {SpWWext,SpWVext,PWext}; 
    FX = SpXext;

    case 'TMat'
    % f(r,n) exterior (Traction)
    %TV
    TVext = @(n) 2*(-n-2).*(n./((2*n+1).*(2*n+3))); %-2(n+2)*SVeg(n);
    %TW 
    TWext = @(n) ((1+2*n.^2)./(1-4*n.^2)); 
    %TX
    TXext = @(n) -(n+2)./(2*n+1); %-(n+2)*SXeg(n); 

    FV = TVext; 
    FW = TWext; 
    FX = TXext; 
    case 'DMat'
    % f(r,n) exterior (double layer)
    %DV
    DVext = @(n) ((2*n.^2+4*n+3)./((2*n+1).*(2*n+3)));
    %DW
    DWext = @(n) ((2*n.^2-2)./(4*n.^2-1)); 
    %DX
    DXext = @(n) ((n-1)./(2*n+1)); 

    FV = DVext; 
    FW = DWext; 
    FX = DXext;  
    case 'DpMat'
    % f(r,n) exterior (double layer)
    %DV
    DVext = @(n) (-n-2).*((2*n.^2+4*n+3)./((2*n+1).*(2*n+3)));
    %DW
    DWVext = @(n) ((-2*n.*(n-1))./(2*n+1)); 
    DWWext = @(n) -n.*((2*n.^2-2)./(4*n.^2-1)); 
      
    %DX
    DXext = @(n) (-n-1).*((n-1)./(2*n+1)); 
    
    %Pressure 
    PWext = @(n) -n.*(n+1); 

    FV = DVext; 
    FW = {DWWext,DWVext,PWext}; 
    FX = DXext;  
    case 'SDMat'
    % f(r,n) exterior (single + double layer)
    %SDV
    SDVext = @(n) (n+1)./(2*n+1);
    %SDW
    SDWext = @(n) (n+1)./(2*n+1); 
      
    %SDX
    SDXext = @(n) n./(2*n+1); 

    FV = SDVext; 
    FW = SDWext; 
    FX = SDXext;  
end
else

switch type 
    case 'SMat'
    % f(r,n) exterior (SMat) 
    FV = @(n) n./((2*n+1).*(2*n+3)); %SVeg; 
    FW = @(n) (n+1)./((2*n+1).*(2*n-1)); %SWeg; 
    FX = @(n) 1./(2*n+1); %SXeg; 
    case 'SpMat'
    % f'(r,n) interior (SpMat)
    SpVVint = @(n) (n+1).*(n./((2*n+1).*(2*n+3))); %)SVeg(n).*r.^n; 
    SpVWint = @(n) (2*(n+1)./(4*n+2));
    SpWint  = @(n) (n-1).*((n+1)./((2*n+1).*(2*n-1))); %(n-1)*SWeg(n).*r.^(n-2);
    SpXint  = @(n) n./(2*n+1); %n*SXeg(n);
    PVint = @(n) -(n+1);

    FV = {SpVVint,SpVWint,PVint}; 
    FW = SpWint; 
    FX = SpXint;

    case 'TSMat'
    % f(r,n) interior (Traction of SL)
    %TV
    %Simplified these two as f(n)r^n+g(n)r^(n-2)  
    TVint = @(n) ((3+4*n+2*n.^2)./(3+8*n+4*n.^2)); 
    %TW
    TWint = @(n) 2*(n-1).*((n+1)./((2*n+1).*(2*n-1))); %SWeg(n); 
    %TX
    TXint = @(n) (n-1)./(2*n+1); %*SXeg(n);

    FV = TVint; 
    FW = TWint; 
    FX = TXint;
    case 'DMat'
    % f(r,n) interior (Double Layer)
    %DV 
    DVint = @(n) ((-2*n.*(n+2))./((2*n+1).*(2*n+3))); 
    %TW
    DWint = @(n) -((2*n.^2+1)./((2*n+1).*(2*n-1))); 
    %TX
    DXint = @(n) -((n+2)./(2*n+1));

    FV = DVint; 
    FW = DWint; 
    FX = DXint;    
    case 'DpMat'
    % f(r,n) interior (normal derivative Double Layer)
    %DV 
    DVVint = @(n) (n+1).*((-2*n.*(n+2))./((2*n+1).*(2*n+3))); 
    DVWint = @(n) -((2*(n+1).*(n+2))./(2*n+1)); 
    %TW
    DWint = @(n) -(n-1).*((2*n.^2+1)./((2*n+1).*(2*n-1))); 
    %TX
    DXint = @(n) -n.*((n+2)./(2*n+1));
    
    PVint = @(n) n.*(n+1);

    FV = {DVVint,DVWint,PVint}; 
    FW = DWint; 
    FX = DXint; 
    case 'SDMat'
    % f(r,n) interior (SL + DL)
    %[S+D]V 
    SDVint = @(n) -n./(2*n+1);  
    %[S+D]W
    SDWint = @(n) -n./(2*n+1); 
    %[S+D]X
    SDXint = @(n) -(n+1)./(2*n+1);

    FV = SDVint; 
    FW = SDWint; 
    FX = SDXint;
end
end

end
