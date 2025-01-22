function [pars,ii2,total_out]=call_randspheroids(N,p,dbox)
    tol=1e-10; comp=false;

    [x1,x2,dst,nghlist,pars,cts,dst2,near,far] = in_sph_MovingBallAlgo(N,tol,p,comp,dbox);
    fprintf("\n Moving Ball Algorithm test for %d spheroids done \n",N);

    A = (dst>0);
    m = min(100,N); 
    m = min(N,size(A,1));
    if m==0
        fprintf("\n no viable spheroids inside constraint");
        return;
    end
    try
        MC1 = maximalCliques(A(1:m,1:m));
        fprintf("\n Max cliques for %d spheroids done \n",m);
    catch ME
        fprintf("\n MaxCliques timed out.")
        error("timed out");
    end

    
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

    total_out=length(ii2);

end