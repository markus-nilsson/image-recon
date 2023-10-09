classdef op_obj_kernel < op_obj % base class for e.g. dictionary actions

    properties
        M;
    end
    
    methods
        
        % Forward model is
        %
        % Y_il  = T_ij w_jk T_kl
        % 
        % where T_ij is the image sampling and T_kl the model kernel

        function O = op_obj_kernel(M)
            if (nargin == 0), return; end
            O.M = M;
            O.n_l = size(M,2); % num experimental contrasts
            O.n_k = size(M,1); % num model coefficients
        end   

        function x = init_x(O, a, b)
            error('not implemented');
        end

        function y = apply(O, x, ind)

            if (nargin < 3), ind = []; end

            if (isnumeric(x))
                y = x * O.M;
            elseif (my_isa(x, 'do_w'))
                y = x.new(x.w * O.M);
            elseif (my_isa(x, 'do_c')) 
                error('not defined'); % need all high res in same space/dim
            else
                error('not defined');
            end

        end

        function y = apply_adjoint(O,x)

            if (isnumeric(x))
                y = x * O.M';
            elseif (my_isa(x, 'do_w'))
                y = x.new(x.w * O.M');
            elseif (my_isa(x, 'do_c'))
                error('not defined');
            else
                error('not defined');
            end
        end

    end

end
