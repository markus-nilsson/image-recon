classdef mec_points_all < mec_points

    methods

        function o = mec_points_all(data)
            o = o@mec_points(data);
            o.coords = o.ind2coords(o.get_ind(1));
        end

        function ind = get_ind(o, ~)
            ind = (1:prod(o.sz))';
        end

        function data_out = apply(o, data, params, ~)

            data_out = data.copy();

            for c_vol = 1:data.n_vol
                data_out.w(:,c_vol) = o.apply_one_vol(data, c_vol, params);
            end

        end    

        function out = calc_displacements(o, params)

            out = o.data.copy();

            for c_vol = 1:o.data.n_vol

                c = o.distort(params, c_vol);
                
                dx = c.x - o.coords.x - o.sz(1)/2;
                dy = c.y - o.coords.y - o.sz(2)/2;
                dz = c.z - o.coords.z - o.sz(3)/2;

                out.w(:,c_vol) = sqrt(dx(:).^2 + dy(:).^2 + dz(:).^2);

            end

        end             

        function out = calc_total_displacements(o, params)

            out = o.calc_displacements(params);
            out.new(sqrt(mean(out.w.^2, 2)));

        end        

    end

end