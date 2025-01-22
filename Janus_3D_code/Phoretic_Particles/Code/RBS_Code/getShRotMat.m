function Rall = getShRotMat(p)

  nmat = floor(p/2)+1;
  Rall = cell(nmat,1);
  matSize =@(n) n*(4*n^2-1)/3;
  
  pstr = ceil(abs(log10(nmat)));
  pstr = ['%0' num2str(pstr) 'd'];
  for nt=1:nmat 
    R = [];
    R = cell(p+1,1);
    fileName = ['rotMatSh_' num2str(p) '_'];
    count = matSize(p+1); 

    RM = readData(fileName, num2str(nt, pstr), 2*count);
    if(~isempty(RM))
      rReal = RM(1:count);
      rImag = RM(count+1:end);
      for n=0:p
        ind = matSize(n);
        dim = (2*n+1);
        R{n+1} = reshape(complex(rReal(ind+1:ind+dim^2), ...
                                 rImag(ind+1:ind+dim^2)),dim,dim);
      end
    else
      u = gl_grid(p);
      theta = u(nt);
      printMsg('* Generating the fast rotation data for p=%d and target=%d.\n',p,nt);

      shcIn = zeros((p+1)*(2*p+1),2*p+1);
      for m=1:2*p+1
        ind = (m-1)*(p+1)+1;
        shcIn(ind:ind+p,m) = 1;
      end
            
      fIn = shSyn(shrinkShVec(shcIn));
      [~, shcOut]= movePole(fIn,theta,0);
      
      for n=0:p
        row = n^2+1:(n+1)^2;
        R{n+1} = shcOut(row,p+1-n:p+1+n);
      end
      
      rReal =[]; rImag = [];
      for n=0:p
        rReal = [rReal;real(R{n+1}(:))];
        rImag = [rImag;imag(R{n+1}(:))];
      end
      writeData(fileName,num2str(nt, pstr),[rReal;rImag]);
    end
    Rall{nt} = R;
  end
  printMsg('* The fast rotation data for p=%d were generate/read from file.\n',p);
