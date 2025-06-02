classdef cost_sense < cost_admm

    properties
        z = 0;
    end

    methods

        function C = cost_sense(mu, ind)
            if (nargin < 2), ind = []; end
            C = C@cost_admm(mu, ind);
        end

        function [f,b] = do_iter(C, x)


            % flat imaginary image
            f = @(x) C.mu * 1i * x.new(imag(x.w));
            b = 0;
   

        end
    end
end