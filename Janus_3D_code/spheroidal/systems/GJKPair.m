function [x1 x2 d] = GJKPair(par1, par2, tol, maxIter)

C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

%check if row vectors
if isequal(size(C1),[1 3])
    C1 = C1.';
    C2 = C2.';
end

iter = 0;
%v0 = [2*x0 2*b1 2*b1].' this is for testing;
v0 = C1 - C2;
v0Norm = norm(v0);
p0 = ellipsoidSupportMapping(-v0, a1, b1, c1, C1, R1);
q0 = ellipsoidSupportMapping(v0, a2, b2, c2, C2, R2);
w0 = p0 - q0;
mu0 = 0;
indices = dec2bin(0:2^4 - 1) - '0';
DeltaX = zeros(16,4);
DeltaX(2,4) = 1;
DeltaX(3,3) = 1;
DeltaX(5,2) = 1;
DeltaX(9,1) = 1;

while (v0Norm - mu0) > tol && iter < 3
    if iter == 0
        v1 = w0;
        pp1 = p0;
        qq1 = q0;
        W1 = w0;
        P1 = p0;
        Q1 = q0;
    else 
        P = [P1 p1];
        Q = [Q1 q1];
        Y = [W1 w1];
        
        cardY = size(Y,2);
        for i=1:cardY - 1
            for j = 1:16
                inZ = find(indices(j,:));
                if sum(indices(j,:)) == i && max(inZ) <= cardY
                    condition1 = true;
                    condition2 = true;
                    k = find(indices(j,:), 1, 'first');
                    notInZ = find(~indices(j,1:cardY));
                    newCondition2 = true;
                    for s = 1:length(notInZ)
                        sumLoop = 0;
                        for t = 1:length(inZ)
                            sumLoop = sumLoop + DeltaX(j,inZ(t))*(dot(Y(:,k), Y(:,inZ(t))) - dot(Y(:,notInZ(s)), Y(:,inZ(t))));

                        end
                        disp(sumLoop);
                        disp(s);
                        newIndex = notInZ(s);
                        newIndices = indices(j,:);
                        newIndices(newIndex) = 1;
                        DeltaX(bin2dec(num2str(newIndices)) + 1,notInZ(s)) = sumLoop;
                        if sumLoop <= 0
                            newCondition2 = true;
                        else
                            newCondition2 = false;
                        end
                        condition2 = newCondition2 && condition2;
                    end 
                    newCondition1 = true;
                    for t = 1:length(inZ)
                        if DeltaX(j,inZ(t)) > 0
                            newCondition1 = true;
                        else
                            newCondition1 = false;
                        end
                        condition1 = newCondition1 && condition1;
                    end
                else
                    continue;
                end

                if condition1 && condition2 
                    v1 = zeros(3,1);
                    for t = 1:length(inZ)
                        v1 = v1 + DeltaX(j, inZ(t)).*Y(:, inZ(t));
                    end
                    v1 = v1./sum(DeltaX(j, :));
                    disp(v1);
                    W1 = Y(:, inZ);
                    P1 = P(:, inZ);
                    Q1 = Q(:, inZ);
                    break;
                end
            
            if condition1 && condition2 
                break;
            end
            end
        end
    end

    v1Norm = norm(v1);
    if v1Norm == 0
        x1 = zeros(3,1);
        x2 = zeros(3,1);
        for t=length(inZ)
            x1 = x1 + DeltaX(j, inZ(t)).*P(:, inZ(t));
            x2 = x2 + DeltaX(j, inZ(t)).*Q(:, inZ(t));
        end
        x1 = x1./sum(DeltaX(j, :));
        x2 = x2./sum(DeltaX(j,:));
        d = 0; 
        return
    end
    p1 = ellipsoidSupportMapping(-v1, a1, b1, c1, C1, R1);
    q1 = ellipsoidSupportMapping(v1, a2, b2, c2, C2, R2);
    w1 = p1 - q1;
    delta1 = dot(v1, w1)/v1Norm;
    mu1 = max([mu0, delta1]);
    iter = iter + 1;
    W0 = W1;
    v0Norm = v1Norm;
    w0 = w1;
    mu0 = mu1;

end

x1 = zeros(3,1);
x2 = zeros(3,1);
for t= 1:length(inZ)
    x1 = x1 + DeltaX(j, inZ(t)).*P(:, inZ(t));
    x2 = x2 + DeltaX(j, inZ(t)).*Q(:, inZ(t));
end 
x1 = x1./sum(DeltaX(j, :));
x2 = x2./sum(DeltaX(j,:));
d = v0Norm;

end
