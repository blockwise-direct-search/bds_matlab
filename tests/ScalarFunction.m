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
                % We set a upper bound for the seed to avoid overflow. We get
                % this upper bound by testing nomad (https://github.com/bbopt/nomad).
                seed = min(abs(ceil(1e5*sin(1e9*sum(x)))) + ...
                       abs(ceil(1e4 * sin(1e7*k_run))) + 5000 * k_run, 2^32 - 1);
                rng(seed)
                if strcmpi(options.noise_type, "uniform")
                    noise = rand(1);
                elseif strcmpi(options.noise_type, "gaussian")
                    noise = randn(1);
                else
                    error("Unknown noise type")
                end
                if options.is_abs_noise
                    f = f + options.noise_level*noise;
                else
                    f = f * (1.0 + max(abs(f), 1) * options.noise_level * noise);
                end
                % If with_gradient is true, it means that we are calculating the fhist of 
                % fminunc and the problem is noisy. 
                % In this case, we provide fminunc with an approximate gradient obtained by
                % finite difference, the step size for the finite difference being 
                % h = sqrt(eps_f), where eps_f is an estimation of the noise in f. 
                % This value of h will improve the performance of fminunc, which 
                % approximates the gradient by finite difference with step size 
                % h_default = sign(x).*sqrt(eps*|x|) when no gradient is provided. 
                % Indeed, the step size that minimizes the error of the finite difference is 
                % h_optimal = sqrt(eps_f/|f''|), where f'' is the second derivative of f, 
                % and h is an approximation of h_optimal when no second-order information
                % is available.
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
                        % Record the function evaluations used by the finite difference in the fhist. 
                        if obj.storeHist
                            obj.valHist(end+1) = f_fd;
                        end
                        obj.nEval = obj.nEval+1;
                        g(i) = (f_fd - f)/h;
                    end

                    % fminunc does not accept the gradient that contains NaN or Inf. Remove such values.
                    g(isnan(g)) = 0;
                    grad_max = 10^10;
                    g = min(grad_max, max(-grad_max, g)); 
                end
                % % Why we need to set f to a huge value when it is NaN or Inf?
                % % The reason is that the algorithm may crash if f is NaN or Inf.
                % if isfield(options, "solver") && strcmpi(options.solver, "nomad")
                %     if isnan(f)
                %         f = min([f, 10^30, sqrt(realmax())]);
                %     end
                % end
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
