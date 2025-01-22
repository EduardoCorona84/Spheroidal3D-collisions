function PF = ShProj(F,isReal)

if nargin==1
    isReal=false; 
end

PF = shSyn(shAna(F),isReal); 

end