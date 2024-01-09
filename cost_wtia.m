classdef cost_wtia < cost_admm

    properties
        z = 0;
    end

    methods

        function C = cost_wtia(mu, ind)
            if (nargin < 2), ind = []; end
            C = C@cost_admm(mu, ind);
        end

        function [f,b] = do_iter(C, x)

            % desired value is zero with this penalty
            tmp = x.w(:, C.ind);

            [max_val,max_ind] = max(tmp, [], 2);

            all_ind = repmat(1:size(tmp,2), size(tmp,1), 1);

            tmp(all_ind ~= max_ind) = 0;

            tmp2 = x.w;
            tmp2(:, C.ind) = tmp;

            C.z = x.new(tmp2);

            [f,b] = C.admm_update(x, C.z);

        end
    end
end