
function y = Lapp(A,x)
if isnumeric(A)
    y=A*x; 
else
    y=real(A(x)); 
end
end
