classdef mec_data < do_w_image_vector

    properties
        interpol_method = 'linear';
        interpolants;
    end

    methods

        function o = mec_data(data, interpol_method)

            o = o@do_w_image_vector(data.w, data.h, data.xps);

            if (nargin > 1)
                o.interpol_method = interpol_method;
            end

            % prepare interpolation kernels
            for c_vol = 1:o.n_vol
                o.interpolants{c_vol} = griddedInterpolant(...
                    o.imreshape(c_vol), ...
                    o.interpol_method, 'nearest');
            end
        end

        function values = interpolate(o, x, y, z, c_vol)

            values = o.interpolants{c_vol}(x, y, z);

        end

        function mask = mask_for_mec(o)

            % Mask the data
            mask = imop_mask(o).imreshape();

            % Expand mask - edges are important
            mask = mio_mask_expand(mask, 5);

            % Trim edges
            mask = o.trim_mask(mask);

        end

        function mask = mask_from_r(o, data_corr, frac)

            % Compute residuals from moving data (see where action is)
            r = o.new(o.w - data_corr.w);
            r = r.new(sqrt(mean( ( ((data_corr.w>0.01) & (o.w>0.01)) .* r.w).^2,2)));
            
            r = mio_smooth_4d(r.imreshape(), 1.1);
            

            % select the most important residual points
            mask = r > quantile(r(:), 1-frac);

            mask = o.trim_mask(mask);

        end

        function mask = trim_mask(o, mask)

            % Trim edges
            mask(1:4, :, :) = 0;
            mask(:, 1:4, :) = 0;
            mask(:, :, 1) = 0;

            mask((end-3):(end), :, :) = 0;
            mask(:,(end-3):(end), :) = 0;
            mask(:,:, end) = 0;
            
        end

    end

end