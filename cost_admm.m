classdef cost_admm < handle 

    properties

        mu1;
        u;
        u_old;

        is_initialized = 0;
        do_mu_update = 0;

    end

    methods

        function C = cost_admm(mu1)
            
            if (nargin < 1), mu1 = 0.1; end            
            C.mu1 = mu1;

        end


        function [f,b] = admm_update(C,x,z)

            % first time initialization
            if (~C.is_initialized)
                C.u = x.zeros();
                C.is_initialized = 1;
            end

            C.u_old = C.u;

            C.u = C.u + (x - z);

            
            b = C.mu1 * (z - C.u);
            f = @(x) C.mu1 * x;

            C.mu_update(x); 
            
        end

        
        function C = mu_update(C, x)
            % Dynamic update of mu according to Boyd et al. 2011

            if (~C.do_mu_update), return; end

            mu = C.mu1;
            y = C.u;
            y_old = C.u_old;

            rs = norm(x - y);
            ss = norm(mu * (y - y_old));
            if rs > 10 * ss
                mu = 2*mu;
                y = y/2;
            elseif ss > 10 * rs
                mu = mu/2;
                y = 2*y;
            end

            C.u = y;
            C.mu1 = mu;

        end

    end
end