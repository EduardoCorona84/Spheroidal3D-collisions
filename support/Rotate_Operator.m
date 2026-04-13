function T = Rotate_Operator(T0,Q,np)

J = [1 3 2]; 
T = Q*reshape(permute(reshape(T0,[3*np 3 np]),J),3,[]);
T = reshape(T,[],3)*Q'; 
T = reshape(permute(reshape(T,[3*np np 3]),J),3*np,3*np); 

end