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

            x_tmp = x;
           
            if (~isempty(C.u))
                x_tmp = x_tmp + C.u; % see Eq. 10 little engine
            end            

            C.z = C.image_filter(x_tmp, C.ind);

            [f,b] = C.admm_update(x, C.z);


        end

        function x_flt = image_filter(C, x, ind)

            if (isempty(ind)), ind = (1:x.h.dim(5)) > 0; end

            TMP2 = x.imreshape();

            for c = find(ind)

                TMP = TMP2(:,:,:,c);

                switch (C.c_type)
                    case 1
                        TMP = medfilt3(real(TMP), [3 3 3]);
                    case 2
                        TMP = mio_smooth_4d(TMP, 0.6);
                    case 3
                        TMP = medfilt3(real(TMP), [5 5 5]);
                    case 4
                        TMP = medfilt3(real(TMP), [7 7 7]);

                    case 5

                        f = @(x) x .* (abs(x) > 10);

                        for k = 1:size(TMP, 3)
                            [wc,ws] = wavedec2(real(TMP), 8, 'db4');
                            TMP(:,:,k) = waverec2(f(wc),ws,'db4');
                        end

                    case 6

                        TMP = ...
                            medfilt3(real(TMP), [3 3 3]) + ...
                            1i * medfilt3(imag(TMP), [3 3 3]);

                    case 7
                        
                        if (c == 1) || (c == 17) || (c == 33)
                            TMP(TMP < 0) = 0;
                        end

                        TMP = medfilt3(TMP, [3 3 3]);

                    case 8

                        %TMP = medfilt3(real(TMP), [3 3 3]);
                        %TMP(TMP < 0) = 0;
                        1;

                end

                TMP2(:,:,:,c) = TMP;

            end

            if (C.c_type == 8)

                S0 = real(TMP2(:,:,:,1));
                S0(S0 < 0) = 0;

                MD = abs(log(TMP2(:,:,:,1)) - log(TMP2(:,:,:,2)));
                MD(isnan(MD(:))) = 0;

                MD(MD > 3.5) = 3.5;
                MD(MD < 0) = 0;
                MD(S0 == 0) = 0;

                % like a soft thresholding
                MD = 0.99 * MD;

                % background removal; power filter
                [rx,ry,rz] = meshgrid(1:5); 
                r2 = (rx-3).^2 + (ry-3).^2 + (rz-3).^2;
                alpha = 0.005;
                filter = 1./(r2 + alpha);
                filter = filter / sum(filter(:));
                MD = imfilter(MD, filter);

                S0_tmp = exp(mean(log(TMP2),4)) .* exp(0.5 * MD);
                S0_tmp = real(S0_tmp);

                S0_tmp(S0_tmp == 0) = S0(S0_tmp == 0);

                S0 = 0.9 * S0 + 0.1 * S0_tmp;

                % apply median filtering last
                S0 = medfilt3(S0, [3 3 3]);
                MD = medfilt3(MD, [3 3 3]);

                TMP2 = cat(4, S0, S0 .* exp(-MD));



            end

            x_flt = x.new(TMP2);
            
        end


    end
end