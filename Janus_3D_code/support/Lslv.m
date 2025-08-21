
function x = Lslv(A,b,parslv)
%global prec;
prec = parslv.prec; 

if nargin<6
   if ~isempty(prec) 
       pr=prec;  
   else
       pr=[]; 
   end
end

if isnumeric(A)
    x=A\b; 
else
    x = zeros(size(b)); 
    %b = parslv.Pr(b); 
    %A = @(x) parslv.Pr(A(x)) + 0.5*(x-parslv.Pr(x)); 
    
    for i=1:size(b,2)
    if nargin==2
        [x(:,i),~,rs,it]=gmres(A,b(:,i),1,1e-6,200,pr);
    else
        [x(:,i),~,rs,it]=gmres(A,b(:,i),parslv.rst,parslv.tol,parslv.maxit,pr);
    end
       fprintf('\n gmres %d iters=%d, res=%1.4g \n',i,prod(it),rs); 
    end
end
end
