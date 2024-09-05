classdef op_obj_kernel_dict_dti < op_obj_kernel_dict

    properties
        xps;

    end
    
    methods
        
        function O = op_obj_kernel_dict_dti(xps, di_list, dd_list, u, n_comp)

            if (nargin < 2), di_list = []; end
            if (nargin < 3), dd_list = []; end
            if (nargin < 4), u = []; end
            if (nargin < 5), n_comp = 16; end
            
            
            if (isempty(di_list)), di_list = linspace(0.3, 3.3, 20) * 1e-9; end
            if (isempty(dd_list)), dd_list = linspace(-0.5, 1, 16); end
            if (isempty(u)), u = uvec_elstat(100, 'froeling'); end

            S_k = zeros(xps.n, numel(di_list), numel(dd_list), size(u,1));
            for i = 1:numel(di_list)
                for j = 1:numel(dd_list)
                    for k = 1:size(u,1) 
                        dt = tm_tpars_to_1x6(3*di_list(i), dd_list(j), u(k,:));
                        S_k(:,i,j,k) = exp(-dt * xps.bt');
                    end
                end
            end

            O = O@op_obj_kernel_dict(S_k, {di_list, dd_list, u}, n_comp);
            O.xps = xps;

        end

        function [mp,s_fit,res] = quick_dti(O,data)

            X = cat(2, ones(O.xps.n,1), O.xps.bt*1e-9);

            S = log(data.w);
            S = real(S);

            FLT = sum(S,1);

            ind = (~isinf(FLT)) & (~isnan(FLT));

            TMP = inv(X' * X) * X' * S(:,ind);

            mp = zeros(size(data.w,2), 7);
            mp(ind,:) = TMP';

            s_fit = (exp(mp * X'))';
            
            res = data.w - s_fit;
        end

    end
end
