function Ax = Acounter(x, A, id)
persistent matVecCnt
if isempty(matVecCnt)
    matVecCnt = containers.Map('KeyType','int32','ValueType','double');
end
if ~exist('id','var') || isempty(id)
    id = 0;
end
if ~matVecCnt.isKey(id)
    matVecCnt(id) = 0;
end

if ischar(x)
    if strcmpi(x, 'cnt')
        Ax = matVecCnt(id);
        return
    elseif strcmpi(x, 'reset')
        Ax = matVecCnt(id);
        keys = matVecCnt.keys;
        K = numel(keys);
        for k = 1:K 
            key = keys{k};
            matVecCnt(key) = 0;
        end
        return
    else
        assert(false, ['Option: ' x ' not recognized'])
    end
end
Ax = A*x;
matVecCnt(id) = matVecCnt(id) + 1;
end