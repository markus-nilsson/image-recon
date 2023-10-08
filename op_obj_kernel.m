classdef op_obj_kernel < op_obj % base class for e.g. dictionary actions

    properties
        n_ims;
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
            O.n_ims = size(M,2); % max(k) in image_sampler ?!
            O.n_mp = size(M,1);
        end        

        function y = apply(O,x)
            y = op_obj.f(x) * O.M;
        end

        function y = apply_adjoint(O,x)
            y = op_obj.f(x) * O.M';
        end

    end

end
