% test ptcharge field

ptch=[1;2];
Xptch=[0 0 0;1 1 1];
Y=[5 5 5];
E=PtChargeField(ptch,Xptch,Y);

a=1/4/pi/(5*sqrt(3))^3*5+1/4/pi*2/(4*sqrt(3))^3*4;
display(E);
display([a a a]);

function E=PtChargeField(ptch,Xptch,Y)
    M=length(ptch);
    np=size(Y,1);
    Rptch=zeros(np,M);
    E=zeros(np,3,M);
    for i=1:M
        % r=Xptch(i,:)-Y;
        r=Y-Xptch(i,:);
        Rptch(:,i)=sqrt(sum(r.^2,2));
        E(:,:,i)=ptch(i)/4/pi./(Rptch(:,i).^3).*r;
    end
    E=sum(E,3);
end