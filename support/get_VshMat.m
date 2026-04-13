function T = get_VshMat(p,type)

sp=(p+1)^2;
VshMattype = ['VshMat_' type]; 

%Load attempt
%printMsg('\n * Reading the Vsh conversion matrix for p=%d \n: ',p);
rFlag=1; 
if p>=9
[TSz,rFlag] = readData([VshMattype '_size'],p,1);
end

if(~rFlag || p<9)

if(p>=9)
printMsg('Failed \n');
printMsg('* Generating the Vsh conversion matrix for p=%d \n',p);
end

inm = @(n,m) n.^2+n+m+1; 
cv = @(m) 1i*m; 
cp1 = @(n,m) (n.*sqrt((n+1).^2-m.^2))./(sqrt((2*n+1).*(2*n+3)));
cm1 = @(n,m) -((n+1).*sqrt(n.^2-m.^2))./(sqrt((2*n+1).*(2*n-1))); 

ii = (1:sp)'; 
nn = floor(sqrt(ii-1)); 
mm=ii-nn.^2-nn-1;

if strcmp(type,'Syn')
%Form T matrix if not available

IGV = ii; JGV=ii+sp; VGV=cv(mm); 
idm = nn>0 & nn>abs(mm);
IGU = ii(idm); JGU=ii(inm(nn(idm)-1,mm(idm))); 
VGU=cm1(nn(idm),mm(idm));
IGU = [IGU;(1:p^2)']; 
JGU = [JGU;inm(nn(1:p^2)+1,mm(1:p^2))]; 
VGU=[VGU;cp1(nn(1:p^2),mm(1:p^2))];

IXU = ii+sp; JXU=ii; VXU=-VGV;
IXV=IGU+sp; JXV=JGU+sp; VXV=VGU; 

JJ=[IGU;IGV;IXU;IXV];
II=[JGU;JGV;JXU;JXV]; 
VV=[VGU;VGV;VXU;VXV];

elseif strcmp(type,'Ana')
    
dn = (nn.*(nn+1))+(nn==0); 
IGV = ii; JGV=ii+sp; VGV=-cv(mm)./dn; 
idm = nn>0 & nn>abs(mm);
IGU = ii(idm); JGU=ii(inm(nn(idm)-1,mm(idm))); 
VGU=cm1(nn(idm),mm(idm))./dn(idm);
IGU = [IGU;(1:p^2)']; 
JGU = [JGU;inm(nn(1:p^2)+1,mm(1:p^2))]; 
VGU=[VGU;cp1(nn(1:p^2),mm(1:p^2))./dn(1:p^2)];

IXU = ii+sp; JXU=ii; VXU=-VGV;
IXV=IGU+sp; JXV=JGU+sp; VXV=VGU; 

II=[IGU;IGV;IXU;IXV];
JJ=[JGU;JGV;JXU;JXV]; 
VV=[VGU;VGV;VXU;VXV];    
    
end

Tdata = [II JJ real(VV) imag(VV)];

if(p>=9)
writeData(VshMattype,p,Tdata);
writeData([VshMattype '_size'],p, size(Tdata,1));
printMsg('* Stored conversion matrix T for p=%d\n', p);
end

else
    
[Tdata,~] = readData(VshMattype,p,[TSz 4]);
    
end

T = sparse(Tdata(:,1),Tdata(:,2),complex(Tdata(:,3),Tdata(:,4)),2*sp,2*sp); 

end