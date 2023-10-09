classdef cost_imfilter_3d < cost_admm

    properties

        c_type;
        z;

    end

    methods

        function C = cost_imfilter_3d(c_type, mu, ind)
            if (nargin < 3), ind = []; end
            C = C@cost_admm(mu, ind);
            C.c_type = c_type;
        end

        function [f,b] = do_iter(C, x)

            C.z = C.image_filter(x, C.ind);

            [f,b] = C.admm_update(x, C.z);

        end

        function x_flt = image_filter(C, x, ind)

            if (isempty(ind)), ind = (1:x.h.dim(5)) > 0; end

            TMP2 = x.imreshape();

            for c = find(ind)

                TMP = TMP2(:,:,:,c);

                switch (C.c_type)
                    case 1
                        TMP = medfilt3(TMP, [3 3 3]);
                    case 2
                        TMP = mio_smooth_4d(TMP, 0.6);
                    case 3
                        TMP = medfilt3(TMP, [5 5 5]);
                end

                TMP2(:,:,:,c) = TMP;

            end

            x_flt = x.new(reshape(TMP2(:), size(x.w)));
        end


    end
end