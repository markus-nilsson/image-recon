classdef cost_imfilter < cost_admm

    properties

        c_type;
        z;
        sz;

    end

    methods

        function C = cost_imfilter(c_type, mu, do_mu_update)
            C = C@cost_admm(mu);
            C.c_type = c_type;    

            if (nargin > 3), C.do_mu_update = do_mu_update; end

        end

        function [f,b] = do_iter(C, x)
            C.z = C.image_filter(x);
            [f,b] = C.admm_update(x, C.z);
        end

        function x_flt = image_filter(C, x)

            assert(isa(x, 'data_obj_image_vector'), 'not an image vector');

            I = x.imreshape();
            for c = 1:size(I,4)
                I(:,:,:,c) = C.apply_filter(I(:,:,:,c));
            end

            x_flt = data_obj_image_volume(I, x.h);

        end

        function x_flt = apply_filter(C, x)            

            switch (C.c_type)
                case 1                    
                    flt = fspecial('gaussian', [5 5], 0.8);
                    x_flt = imfilter(x, flt);
                case 2
                    x_flt = medfilt3(x, [3 3 3]);
                case 3
                    x_flt = TVL1denoise(x, 1); 
                case 4
                    x_flt = medfilt2(real(x), [3 3]) + ...
                        1i * medfilt2(imag(x), [3 3]);
                case 5
                    x_flt = wdenoise(x);


            end

        end



    end
end