classdef cost_positivity < cost_admm

    properties
        z = 0;
    end

    methods

        function C = cost_positivity(mu, ind)
            if (nargin < 2), ind = []; end
            
            C = C@cost_admm(mu, ind);

        end

        function [f,b] = do_iter(C, x)

            % desired value is zero with this penalty
            tmp = x.w;
            tmp(tmp < 0) = 0;
            C.z = x.new(tmp);

            [f,b] = C.admm_update(x, C.z);

        end
    end
end