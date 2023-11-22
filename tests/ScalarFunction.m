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
            %   STOREHIST to decide whether to store the history.
            obj.userFun = varargin{1}.objective;
            obj.storeHist = true;
            if nargin > 1
                obj.storeHist = logical(varargin{2});
            end
            obj.nEval = 0;
            obj.valHist = [];
        end

        function [f, g] = fun(obj,x,is_noisy,k_run,options)
            %FUN evaluates the scalar function. If and only if the
            %   function is noisy and the with_gradient is true, it also
            %   returns the gradient of the function.
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
                % We set nan values to zero since we only need x to define
                % the seed and we are not evaluating any function at x.
                x(isnan(x)) = 0;
                seed = abs(ceil(1e5*sin(1e9*sum(x)))) + ...
                       abs(ceil(1e4 * sin(1e7*k_run))) + 5000 * k_run;
                rng(seed)
                if strcmpi(options.noise_type, "uniform")
                    noise = rand(1);
                elseif strcmpi(options.noise_type, "gaussian")
                    noise = randn(1);
                else
                    error("Unknown noise type")
                end
                if options.is_abs_noise
                    f = f+options.noise_level*noise;
                else
                    f = f*(1.0+options.noise_level*noise);
                end
                % If with_gradient is true, it means that we are calculating the fhist of 
                % fminunc and the problem is noisy. In this case, we offer the gradient
                % by finite difference using another FiniteDifferenceStepSize, which
                % is sqrt(max(abs(f), 1)*options.noise_level). It will improve 
                % the performance of fminunc than using the default FiniteDifferenceStepSize.
                % We need to make sure the function evaluations during the finite difference
                % are included in the fhist. Since fminunc does not accept the gradient that
                % contains nan, we set nan in the gradient to be zero. Since fminunc also does
                % not accept the gradient which contains inf or -inf, we set the gradient to 
                % be 10^10 if it is larger than 10^10 and -10^10 if it is smaller than -10^10 piecewisely. 
                if isfield(options, "with_gradient") && options.with_gradient && nargout >= 2
                    if options.is_abs_noise
                        h = sqrt(options.noise_level);
                    else
                        h = sqrt(max(abs(f), 1)*options.noise_level); 
                    end
                    dim = length(x);
                    g = NaN(dim, 1);
                    V = eye(dim);
                    for i = 1:dim
                        f_fd = obj.userFun(x + h*V(:, i));
                        if obj.storeHist
                            obj.valHist(end+1) = f_fd;
                        end
                        obj.nEval = obj.nEval+1;
                        g(i) = (f_fd - f)/h;
                    end
                    g(isnan(g)) = 0;
                    grad_max = 10^10;
                    g = min(grad_max, max(-grad_max, g)); 
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
