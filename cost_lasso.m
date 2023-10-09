classdef cost_lasso < cost_admm

    properties
        z = 0;
        lambda;
    end

    methods

        function C = cost_lasso(mu, ind, lambda)
            if (nargin < 2), ind = []; end
            C = C@cost_admm(mu, ind);
            C.lambda = lambda;
        end

        function [f,b] = do_iter(C, x)

            if (isempty(C.u))
                tmp = x;
            else
                tmp = x + C.u;
            end
            
            C.z = cost_lasso.soft_thresholding(...
                tmp, C.lambda / C.mu);

            [f,b] = C.admm_update(x, C.z);

        end
    end

    methods (Static)

        function y = soft_thresholding(x, lambda)

            y = x.new(x.w);

            ax = abs(x.w);
            ind = ax < lambda;
                
            y.w = (ax - lambda) .* sign(x.w);
            y.w(ind) = 0;

        end
    end
end