function [ShTg,ShTr,ShY] = get_Surfgrad_coeffs(px,p,sprs)

if sprs
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
       TrData = [TrData ; [ig ; ir ; iy] [jg ; jr ; jy] [real(vg) ; real(vr) ; real(vy)] [imag(vg) ; imag(vr) ; imag(vy)] [ng ; nr ; ny]];  
    end
    
    writeData('TractionCoeffsSp',p, TrData);
    writeData('TractionCoeffsSpSize',p, size(TrData,1));
    printMsg('* Stored generated Traction coefficient matrices for p=%d\n', p);
  else
    sp = (p+1)^2; 
    [TrData,rFlag] = readData('TractionCoeffsSp',p,[TrSz 5]);  
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
  
else
    printMsg('* Reading the Traction kernel coefficient matrices for p=%d: ',p);
    sp = (p+1)^2; 
  [TrData,rFlag] = readData('TractionCoeffsDn',p,[9*sp 18*sp]);

  if(~rFlag)
    printMsg('Failed\n');
    printMsg('* Generating the Traction kernel coefficient matrices for p=%d\n',p);

    [ShTg,ShTr,ShY] = Surfgrad_coeffs(px,p,1,'VW','VW',0);
    
    TrData = [full(real(cell2mat(ShTg))) full(imag(cell2mat(ShTg))) ...
        full(real(cell2mat(ShTr))) full(imag(cell2mat(ShTr))) ...
        full(real(cell2mat(ShY))) full(imag(cell2mat(ShY)))]; 
    
    writeData('TractionCoeffsDn',p, TrData);
    printMsg('* Stored generated Traction coefficient matrices for p=%d\n', p);
  else 
    ShTg = cell(3,1); ShTr = ShTg; ShY = ShTg;  
    
    for j=1:3
        ind = (1+3*sp*(j-1)):3*sp*j;  
        ShTg{j} = complex(TrData(ind,1:3*sp),TrData(ind,3*sp+1:6*sp));   
        ShTr{j} = complex(TrData(ind,6*sp+1:9*sp),TrData(ind,9*sp+1:12*sp));
        ShY{j} = complex(TrData(ind,12*sp+1:15*sp),TrData(ind,15*sp+1:18*sp));
    end
    
    printMsg('Successful\n');
  end
end
end