classdef cost_penalty < cost_admm

    properties
        z = 0;
    end

    methods

        function C = cost_penalty(mu, ind)
            if (nargin < 2), ind = []; end
            C = C@cost_admm(mu, ind);
        end

        function [f,b] = do_iter(C, x)

            % desired value is zero with this penalty
            C.z = x.new(zeros(size(x.w)));

            [f,b] = C.admm_update(x, C.z);

        end
    end
end