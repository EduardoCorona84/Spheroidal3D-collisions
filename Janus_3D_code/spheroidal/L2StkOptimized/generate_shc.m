function varargout = generate_shc(varargin)
    if nargin == 0, error('Inputs must be passed in.'); end
    for i=1:nargin
        varargout{i} = shAna(varargin{i});
    end
end