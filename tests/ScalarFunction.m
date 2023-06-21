classdef ScalarFunction < handle
    %SCALARFUNCTION Scalar function.
    %   This class defines a scalar function f : R^n -> R and stores
    %   information related to its evaluations.


    properties (Access = private)
        userFun  % the function provided by the user
        storeHist  % whether to store the history
    end
    properties
        nEval      % number of function evaluations
        valHist    % history of function evaluations
    end

    methods
        function obj = ScalarFunction(varargin)
            %SCALARFUNCTION Construct an instance of this class.
            %
            %   OBJ = SCALARFUNCTION(P) returns an instance of this class
            %   which stores the history.
            %
            %   OBJ = SCALARFUNCTION(P,STOREHIST) uses the value of
            %   STOREHIST to decide whethe to store the history.
            obj.userFun = varargin{1}.objective;
            obj.storeHist = true;
            if nargin > 1
                obj.storeHist = logical(varargin{2});
            end
            obj.nEval = 0;
            obj.valHist = [];
        end

        function f = fun(obj,x,is_noisy,k_run,options)
            %FUN Evaluation the scalar function.
            %
            %   F = OBJ.FUN(X) returns the function evaluation at X.
            f = obj.userFun(x);
            obj.nEval = obj.nEval+1;
            if obj.storeHist
                obj.valHist(end+1) = f;
            end
            if nargin <= 2
                is_noisy = false;
            end
            if is_noisy
                x(isnan(x)) = 0;
                seed = abs(ceil(1e5*sin(sum(x)))) + ...
                       abs(ceil(1e4 * sin(1e3*k_run))) + 5000 * k_run;
                rng(seed)
                if strcmpi(options.noise_type, 'uniform')
                   noise = rand(1);
                end
                if strcmpi(options.noise_type, 'gaussian')
                   noise = randn(1);
                end
                if options.is_abs_noise
                    f = f+options.noise_level*noise;
                else
                    f = f*(1.0+options.noise_level*noise);
                end
            end
        end

        function valHist = get.valHist(obj)
            %GET.VALHIST Get attribute for VALHIST.
            %
            %    VALHIST = OBJ.VALHIST returns a copy of the history of
            %    function evaluations.
            valHist = obj.valHist(:);
        end
    end
end
