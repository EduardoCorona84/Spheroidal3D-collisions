function flagOut  = printMsg(format, varargin)
% PRINTMSG(format,...) - Controls the destination of the output messages.
%
% EXAMPLES:
%     printMsg('Some text') -- prints the text
%    
%     printMsg([],'verbose','screen') or printMsg([],'verbose','file','fileName')

%     set the verbose option. To set the verbose option the format string
%     should be empty, otherwise printMsg prints out the input arguments
%     according to the given format. The verbose option is either 'quiet',
%     'screen', or 'file'. If it is file, it should be followed by the file
%     name string.    
%

  persistent verbose fileName
  
  %-- Setting the verbose mode
  if(isempty(format) && any(strcmp(varargin,'verbose')))
    ind = find(strcmp(varargin,'verbose'));
    verbose = varargin{ind+1};
    if(strcmp(verbose,'file')), fileName = varargin{ind+2};end 
  end
  if(isempty(verbose)), verbose = 'screen'; end

  %-- separation
  sep = [];
  if(any(strcmp(varargin,'sep')))
    ind = find(strcmp(varargin,'sep'));
    sep = varargin{ind+1};
    varargin = varargin([1:ind-1 ind+2:length(varargin)]);
  end

  flagOut = 0;
  if(isempty(format) && isempty(varargin))
    return
  end

  %-- printing
  switch verbose
   case 'quiet'
    return;
   case 'screen'
    fid = 1;
    flag = 1;
   case 'file'
    fid = fopen(fileName,'a');
    flag = 2;
  end

  str = sprintf(format,varargin{:});
  if(~isempty(sep))
    str = [str '\n' repmat(sep, 1, length(str)) '\n'];
  end
  fprintf(fid,str);
    
  %fprintf(fid,[indent format],varargin{:});
  if(fid~=1), fclose(fid); end
  flagOut = flag;