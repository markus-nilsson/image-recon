classdef cost_admm < handle 

    properties

        mu;
        u;
        u_old;

        is_initialized = 0;
        do_mu_update = 0;

        ind; % 4th indices for filtering
        O; % for selecting indices to add a cost to
        

    end

    methods

        function C = cost_admm(mu, ind)
           
            if (nargin < 1), mu = 0.1; end            
            if (nargin < 2), ind = []; end

            C.mu = mu;
            C.ind = ind;

            if (isempty(C.ind))
                C.O = 1;
            else
                C.O = op_obj_kernel(diag(ind));
            end
            
        end


        function [f,b] = admm_update(C,x,z)

            % first time initialization
            if (~C.is_initialized)
                C.u = x.zeros();
                C.is_initialized = 1;
            end

            % Compared with Assländer's implementation, we have
            % P = 1
            % G = z and C.z
            % z = C.u

            % Compared with Boyd, we have
            % y = u
            % with A = 1, B = -1, c = 0


            C.u_old = C.u;

            C.u = C.u + (x - z);

            
            b = C.mu * (C.O * (z - C.u));
            f = @(x) C.mu * (C.O * x);

            C.mu_update(x); 
            
        end

        
        function C = mu_update(C, x)
            % Dynamic update of mu according to Boyd et al. 2011

            % This does not work at present, don't know why

            if (~C.do_mu_update), return; end


            mu = C.mu;
            y = C.u;
            y_old = C.u_old;

            rs = norm(x - y);
            ss = norm(mu * (y - y_old));
            if rs > 10 * ss
                mu = 2*mu;
                y = y * (1/2);
            elseif ss > 10 * rs
                mu = mu/2;
                y = 2*y;
            end

            

            C.u = y;
            C.mu = mu;

        end

    end
end