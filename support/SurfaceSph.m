classdef SurfaceSph < handle 
        
    properties(SetObservable) %i.e. changes of these will cause other
                              %memebers to get updated

        cart = vec3d;        %Cartesian coordinate of the points. cart is of
                             %type vec3d that is user defined.
        
        upFreqCoeff;     %upsampling rate for differentiation
       
        filterFreqCoeff; %filtering rate for differentiation
        stokesOpStale = true;   %flag indicating whether the precomputed
                                %stokes is stale or not (in case the stokes
                                %operator is stored as a matrix)
        dblLayerOpStale = true; %same as above for double layer
    end
    
    properties(SetObservable, SetAccess = 'protected')
        p             %The number of spherical harmonics frequencies. The
                      %number of points is 2*p*(p+1).
    end
    
    properties(SetAccess = 'protected')
        shc = vec3d;         %Spherical harmonic coefficients of Cartesian
                             %coordinates. shc is of type vec3d (user defined
                             %type).
        
        geoProp;             %Geometric properties of the surface. It
                             %contains first and second derivative, first and
                             %second fundamental coefficient, normal vector,
                             %Gaussian and mean curvature.
        
        upFreq;              % upsampling frequency for differentiation  = upFreqCoeff * p
        filterFreq;          % filtering frequency for differentiation = filterFreqCoeff * p
    end
        
    methods

        %Constructor for the objects of type SurfaceSph. It can be initialized
        %by giving the Cartesian coordinate of the points. Coordinate points
        %are stored in (and can be initialized with) an object of type vec3d.
        function obj=SurfaceSph(cartInit, upFreqCoeffIn, filterFreqCoeffIn)
            
            if(isa(cartInit,'SurfaceSph'))
                obj.upFreqCoeff = cartInit.upFreqCoeff;
                obj.filterFreqCoeff = cartInit.filterFreqCoeff;
                obj.cart = cartInit.cart;
            else               
                if(nargin>1)
                    obj.upFreqCoeff = upFreqCoeffIn;
                else
                    obj.upFreqCoeff = 1;
                end
                if(nargin>2)
                    obj.filterFreqCoeff = filterFreqCoeffIn;
                else
                    obj.filterFreqCoeff = 1;
                end
            
                obj.cart = vec3d(cartInit);
            end
        end
        
        function set.cart(obj,val)
            %Set method for cart field. cart field is of type vec3d and can
            %be set to any object of this type.
            
            if(isa(val,'vec3d'))
                obj.cart = val;
            else
                error(['Surface''s Cartesian coordinates can be set by an object ' ...
                    'of class vec3d.'])
            end
       
            % Updating the rest of the params
            obj.p = (sqrt(2*size(obj.cart.x,1)+1)-1)/2; %#ok<*MCSUP>
            obj.shc = obj.cart.shAna();
            try 
                obj.geoProp = calcGeoProp(obj);
            catch 
                obj.geoProp = [];
            end
            obj.stokesOpStale = true;
            obj.dblLayerOpStale = true;
        end
        
        function set.upFreqCoeff(obj,val)
            if( size( val ) ~= size( obj ) )
              val = repmat( val, size( obj ) );
            end
            for ii=1:length(obj)
              flag = (obj(ii).upFreqCoeff ~= val(ii));
              obj(ii).upFreqCoeff = val(ii); 
              if(flag), obj(ii).resample; end
            end
        end
        
        function set.filterFreqCoeff(obj,val)
            if( size( val ) ~= size( obj ) )
              val = repmat( val, size( obj ) );
            end
            for ii=1:length(obj)
              flag = (obj(ii).filterFreqCoeff ~=val(ii));
              obj(ii).filterFreqCoeff = val(ii);
              if(flag), obj(ii).resample; end
            end
        end
        
        function f = get.upFreq(obj)
            f = ceil(obj.p*obj.upFreqCoeff);
        end
        
        function f = get.filterFreq(obj)
            f = floor(obj.p*obj.filterFreqCoeff);
        end
        
        function [Sf, callsRet] = stokesMatVec(obj, f, varargin)
       % Stokes matvec wrapper     
            persistent calls
            
            if ( isempty(calls) )
                calls.n = 0;
                calls.t = 0;
            end
            
            if(nargin==1)
                Sf = [];
                callsRet = calls;
                return;
            end
            
            tt = clock;
            Sf = kernelS(f, obj, varargin{:});

            calls.t = calls.t + etime(clock,tt);
            if(~isempty(f)), calls.n = calls.n + 1;end
            if(nargout>1), callsRet = calls; end
        end
        
        function [Df,numcalls] = doubleLayerMatVec(obj, f, varargin)
            %Double layer operator wrapper
              
            persistent ncalls;
            if ( isempty(ncalls) ), ncalls = 0; end
            
            ncalls = ncalls + 1;
            numcalls = ncalls;
                       
            Df = kernelD(f, obj, varargin{:});
        end               
        
        function hOut = plot(obj,varargin)
            h=plotb(obj.getAllCart(),varargin{:}); drawnow;
            if(nargout>0), hOut = h; end
        end
        
        function h = quiver(obj,vec,varargin)
            h = quiver3(obj.cart.x,obj.cart.y,obj.cart.z,...
                vec.x,vec.y,vec.z,varargin{:});
        end
            
        function obj = resample(obj,newFreq)
            if(nargin<2), newFreq = obj.p; end
            if(~isempty(obj.cart))
                obj.cart = obj.cart.interpsh(newFreq);
            end
        end
        
        function X = movePole(obj, thetaIdx, phiIdx)
            XX = movePoleShcMatrix(reshape(obj.cart.to_array(), [], 3),...
                thetaIdx, phiIdx);
            X(1:length(XX)) = vec3d;
            for ii=1:length(XX)
                X(ii) = vec3d(XX{ii});
            end
        end
        
        function filter(obj,filtFreq)
            obj.cart = obj.cart.filtersh(filtFreq);
        end
        
        function a = area(obj)
            [~,a] = reducedVolume(obj);
        end
        
        function v = volume(obj)
            v = reducedVolume(obj);
        end
        
        function cent = centerOfMass(obj)
            cent = getCenter(obj);
        end
        
        function ra = reducedVol(obj)
            [~,~,ra] = reducedVolume(obj);
        end
        
        function II = momentOfInertiaVolume(obj)
            C = obj.centerOfMass();
            R = obj.cart-C;
            S = .5*dot(R,R);
            w = .5*R.*R.*R;
            nor = obj.geoProp.nor;
            II = eye(3);
            for i=1:3
                ei=vec3d(II(:,i));                               
                f = integrateOverS(obj, dot(w,nor).*ei - dot(R,ei).*S.*nor);
                II(:,i) = f.to_array;
            end            
            II = II/obj.volume(); %normalize with volume (equivalent to whole mass)
            %expect .4 for sphere and [1 1 .4] for ellipse (1,1,2)
        end
        
        function II = momentOfInertiaShell(obj)
            C = obj.centerOfMass();
            R = obj.cart-C;
            a = dot(R,R);
            II = eye(3);
            for i=1:3
                ei=vec3d(II(:,i)); 
                f = integrateOverS(obj, a*ei-dot(R,ei).*R);
                II(:,i) = f.to_array;
            end            
            II = II/obj.area(); %normalize with area (equivalent to shell mass)
        end
                
        function X = getAllCart(obj)
            for ii=1:length(obj)
                X(ii) = obj(ii).cart;
            end
        end                
    end
end
