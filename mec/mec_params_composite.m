classdef mec_params_composite < mec_params_base

    properties
        p_global
        p_volume
        p_slice
        current
    end

    methods

        function o = mec_params_composite(data, order)
            
            o = o@mec_params_base(data);

            o.p_global = mec_params_global_orders(data, order);
            o.p_volume = mec_params_volume_rigid(data);
            o.p_slice  = mec_params_slice(data);

            o.enabled = true;            

        end

        function o = reset(o)

            if isempty(o.p_global), return; end

            o.p_global = o.p_global.reset();
            o.p_volume = o.p_volume.reset();
            o.p_slice  = o.p_slice.reset();

        end        

        function o = update(o, x, c_vol)

            switch (o.current)
                case 'global'
                    o.p_global = o.p_global.update(x);
                case 'volume'
                    o.p_volume = o.p_volume.update(x, c_vol);
                case 'slice'
                    o.p_slice = o.p_slice.update(x, c_vol);
            end
        end

        function x = get_x(o, c_vol)

            if isempty(o.current)
                x = 1;
                return; 
            end            

            switch (o.current)
                case 'global'
                    x = o.p_global.get_x();
                case 'volume'
                    x = o.p_volume.get_x(c_vol);
                case 'slice'
                    x = o.p_slice.get_x(c_vol);
            end            
        end

        function p = get_p(o, ~)

            switch (o.current)
                case 'global'
                    p = o.p_global.get_p();
                case 'volume'
                    p = o.p_volume.get_p(c_vol);
                case 'slice'
                    p = o.p_slice.get_p(c_vol);
            end            

        end

        function d = distort(o, c, c_vol)

            d = o.p_global.distort(c, c_vol);
            d = o.p_volume.distort(d, c_vol);
            d = o.p_slice.distort(d, c_vol);

        end


        function plot(o)
            1;
        end

    end

    methods

        function o = optimize(o, level, data_mov, data_ref, points)

            disp('--------------------------------------');
            disp(sprintf('Optimizing %s', level));
            disp('--------------------------------------');

            o.current = level;

            switch (o.current)
                case 'global'
                    o.opt_global = true;
                    o.opt_opts.Display = 'iter';
                otherwise 
                    o.opt_global = false;
                    o.opt_opts.Display = 'none';
            end

            o.opt_opts.TypicalX = ones(size(o.get_x(1)));
    
            o = optimize@mec_params_base(o, data_mov, data_ref, points);

        end
    end

end