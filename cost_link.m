classdef cost_link < cost_admm

    properties
        z = 0;
    end

    methods

        function C = cost_link(mu, ind)
            if (nargin < 2), ind = []; end

            if (sum(ind) ~= 2), error('stop'); end
            
            C = C@cost_admm(mu, ind);

        end

        function [f,b] = do_iter(C, x)

            a = min(find(C.ind));
            b = max(find(C.ind));

            sc = x.w(:,a) \ x.w(:,b);

            C.z = x.new(zeros(size(x.w)));
            C.z.w(:,a) = x.w(:,a);
            C.z.w(:,b) = x.w(:,a) * sc;

            [f,b] = C.admm_update(x, C.z);

        end
    end
end