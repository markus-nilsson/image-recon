classdef op_obj_kernel_stack < op_obj

    properties

        K_list;
        n_kernels;
    end
    
    methods
        
        function O = op_obj_kernel_stack(K_list)

            O = O@op_obj();
            O.K_list = K_list;
            O.n_kernels = numel(K_list);
            
        end

        function x = init_x(O)
            x = do_c(O.n_kernels);
            for c = 1:O.n_kernels
                x.data_obj{c} = O.K_list{c}.init_x();
            end
        end

        function y = apply(O, x)

            for c = 1:O.n_kernels
                tmp = O.K_list{c}.apply(x.data_obj{c});
                if (c == 1)
                    y = tmp;
                else
                    y = y + tmp;
                end
            end
        end


        function y = apply_adjoint(O, x)

            y = do_c(O.n_kernels);

            if (x.c_type == 1)
                for c = 1:O.n_kernels
                    y.data_obj{c} = O.K_list{c}.apply_adjoint(x);
                end
            elseif (x.c_type == 2)
                for c = 1:O.n_kernels
                    y.data_obj{c} = O.K_list{c}.apply_adjoint(x.data_obj{c});
                end
            end
            
        end





    end
end
