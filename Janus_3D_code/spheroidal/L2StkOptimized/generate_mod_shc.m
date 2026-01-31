function varargout = generate_shc(Gmtx, varargin)
    if nargin == 1, error('Inputs must be passed in.'); end
    for i=1:nargin-1
        varargout{i} = Gmtx \ varargin{i};
    end
end