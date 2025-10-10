function [val rFlag] = readData(fileInit, suffix, dim)
% READDATA(fileInit, suffix, dim) - Project specific interface to read binary files
try
    root = getenv('SLURM_SCRATCH');
catch
    root = '../data';
end
nd = prod(dim); % number of entries to be read
if(isstr(suffix))
    fileName = [root filesep fileInit suffix '.bin'];
else
    fileName = [root filesep fileInit num2str(suffix) '.bin'];
end
rFlag = exist(fileName,'file');
if(rFlag)
    fid = fopen(fileName,'r');
    val = fread(fid, nd,'double');
    fclose(fid);

    if(size(dim,2)>1)
        val = reshape(val, dim);
    end
else
    val = [];
end
