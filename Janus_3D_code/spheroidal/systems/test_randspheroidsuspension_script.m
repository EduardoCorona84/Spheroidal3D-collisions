 function test_randspheroidsuspension_script()

N = 100; tol=1e-10; p=8; comp=false; 
tic; 
[x1,x2,dst,nghlist,pars,cts,dst2,near,far] = test_MovingBallAlgo(N,tol,p,comp);
fprintf("\n Moving Ball Algorithm test for %d spheroids done \n",N);

pause(0.1); 

A = (dst>0);
m = min(100,N); 
MC1 = maximalCliques(A(1:m,1:m));
fprintf("\n Max cliques for %d spheroids done \n",m);

pause(0.1); 

ii=1:N; 
id1 = MC1(:,1)==1;
ii2 = ii(id1);

if m<N
    for i=m+1:N
    Airow = A(i,ii2); 
    Si = sum(Airow)/length(ii2);
        if abs(Si-1)<1e-8
            ii2 = [ii2 i];
            display(length(ii2));
        end
    end
end

fprintf("\n Ended up with %d spheroids \n",length(ii2));

% Plot suspension
p=16; 
Surf{1} = SurfaceSph(shape_gallery(p,'ellipseZ'));
Surf{2} = SurfaceSph(shape_gallery(p,'oblateZ'));
Xs=cell(2,1); 
for k=1:2
Xs{k} = reshape(Surf{k}.cart.to_array,[],3);
end
np=2*p*(p+1);

for i=1:length(ii2)
    if strcmp(pars(i).shape,'ellipseZ')
        X = Xs{1};
    else
        X = Xs{2};
    end
    Xi = (pars(i).R*X')' + repmat(pars(i).C,np,1);
    plotb(vec3d(Xi)); hold on;
end

% D = dst(ii2,ii2);
% figure; 
% histogram(min(D+10*eye(size(D,1)))); 
% figure; 
% histogram(log10(min(D+10*eye(size(D,1)))));

end