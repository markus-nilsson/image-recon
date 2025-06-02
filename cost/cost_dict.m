classdef cost_dict < cost_admm

    properties
        
        O;
        D;
        z;

        data = [];

    end

    methods

        function C = cost_dict(O, mu)
            C = C@cost_admm(mu);
            C.O = O;
        end
 
        function [f,b] = do_iter(C, x)

            if (all(x(:) == 0))
                f = @(x) 0;
                b = 0;
                return;
            end

            C.D = C.O.map_to_dict(x);
            C.z = C.proj_to_dict(x);

            % should not be necessary
            if (~isempty(C.data))
                rd = C.O * C.z;
                C.z = C.z / ((C.data.w(:)\rd(:)) + eps);
            end

            [f,b] = C.admm_update(x, C.z);

        end

        % Project into dictionary
        function x = proj_to_dict(C, x)
            x = C.D .* (sum(C.D .* x,1));

            x(isnan(x(:))) = 0;
        end


    end
end