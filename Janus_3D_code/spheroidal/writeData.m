function count = writeData(fileName, suffix, data)
% WRITEDATA(fileName, data) - Project specific interface to write binary files
this_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(this_dir, 'data');

if (ischar(suffix))
    fileName = fullfile(data_dir, [fileName suffix '.bin']);
else
    fileName = fullfile(data_dir, [fileName num2str(suffix) '.bin']);
end

if(~exist(data_dir, 'dir'))
    mkdir(data_dir) 
end

fid = fopen(fileName, 'w');
if(fid>0)
	count = fwrite(fid, data, 'double');
	fclose(fid);
else
	count = fid;
end
