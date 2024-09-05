classdef cost_multi_dict < cost_admm

    properties
        
        D;
        z;
        K;

        data = [];

        c_iter;

    end

    methods

        function C = cost_multi_dict(K, mu)

            if (~isa(K, 'op_obj_kernel_stack'))
                error('expected op_obj_kernel_stack');
            end

            C = C@cost_admm(mu);
            C.O = 1;
            C.K = K;
            C.c_iter = 0;
        end
 
        function [f,b] = do_iter(C, x)

            C.c_iter = C.c_iter + 1;

            if (C.c_iter == 1) % this is a bit stupid
                C.z = C.K.init_x();
                [f,b] = C.admm_update(x, C.z);
                return;
            end

            
            for c = 1:numel(C.K.K_list)
                C.D{c} = C.K.K_list{c}.map_to_dict(x.data_obj{c});
                C.z.data_obj{c} = C.proj_to_dict(x.data_obj{c}, C.D{c});
            end

%             C.D = C.O.map_to_dict(x);
%             C.z = C.proj_to_dict(x);
% 
%             % should not be necessary
%             if (~isempty(C.data))
%                 rd = C.O * C.z;
%                 C.z = C.z / ((C.data.w(:)\rd(:)) + eps);
%             end

            [f,b] = C.admm_update(x, C.z);

        end

        % Project into dictionary
        function x = proj_to_dict(C, x, D)
            TMP = D.w .* (sum(D.w .* x.w,2)) ./ sum(D.w .* D.w, 2);

            TMP(isnan(TMP(:))) = 0;

            x = x.new(TMP);
        end


    end
end