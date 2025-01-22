function GradF = Sph_spdiff(F,isReal)

if nargin==1
    isReal=false; 
end

[d1,d2] = size(F);
p = (sqrt(2*d1+1)-1)/2;
sp = (p+1)^2; 

shF = shAna(F); 

zsh = zeros(sp,d2); 
shRGX = [zsh;shF;zsh]; 
GradF = VshSyn(shRGX,'GX',isReal); 

end