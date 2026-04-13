function count = writeData(fileName, suffix, data)
% WRITEDATA(fileName, data) - Project specific interface to write binary files
global DATA_DIR

if ~isempty(DATA_DIR) 
    root = DATA_DIR;
else
    root = '../data';
end
if(isstr(suffix))
  fileName = [root filesep fileName suffix '.bin'];
else
  fileName = [root filesep fileName num2str(suffix) '.bin'];
end

if(~exist(root, 'dir'))
  mkdir(root)
end

fid = fopen(fileName, 'w');
if(fid>0)
	count = fwrite(fid, data, 'double');
	fclose(fid);
else
	count = fid;
end
