function Ax = Acounter(x, A, transpose)
persistent matVecCnt

if isempty(matVecCnt)
    matVecCnt = 0;
end
if ~exist('transpose','var') || isempty(transpose)
    transpose = false;
end

if ischar(x)
    if strcmpi(x, 'cnt')
        Ax = matVecCnt;
        return
    elseif strcmpi(x, 'reset')
        Ax = matVecCnt;
        matVecCnt = 0;
        return
    else
        assert(false, ['Option: ' x ' not recognized'])
    end
end
if transpose
    Ax = A'*x;
else
    Ax = A*x;
end
matVecCnt = matVecCnt + 1;
end