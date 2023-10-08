classdef op_obj_append < op_obj

    properties

        A, B;

    end

    properties (Access = private)
        c_type;
    end
    
    methods
        
        function O = op_obj_append(A,B)
            
            % assume A and B are of types kernel or image sampler

            type_a = 'op_obj_image';
            type_b = 'op_obj_kernel';

            if (0) 
                1;
            elseif (my_isa(A, type_a) && my_isa(B, type_a))
                O.c_type = 1;
            elseif (my_isa(A, type_b) && my_isa(B, type_b))
                O.c_type = 2;
            elseif (my_isa(A, type_a) && my_isa(B, type_b)) 
                O.c_type = 3;
            elseif (my_isa(A, type_b) && my_isa(B, type_a)) 
                O.c_type = 4;
            else 
                error('check this');
            end

            switch (O.c_type)
                case 1 % two image samplers
                    O.n_mp  = max(A.n_mp, B.n_mp);
                    O.n_vox = A.n_vox;
                    O.k = A.k;
                case 2 % two dictionaries – no use case?
                    error('check this');
                case 3 % image sampling first, then addition - probably slower
                    error('check this'); 
                case 4
                    O.n_mp = A.n_mp;
                    O.n_vox = B.n_vox;
            end

            O.A = A;
            O.B = B;

        end

        function x = init_x(O)
            x = O.A.init_x(O.n_vox, O.n_mp);
        end
        
        function y = apply(O,x)
            y = O.B * (O.A * x);
        end

        function y = apply_adjoint(O,x)
            y = O.A' * (O.B' * x);
        end

    end
end
