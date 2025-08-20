classdef mec_points

    properties
        n_vol;
        sz;
        n;
        data;
        c_vols = [];

        ind;
        coords; % x, y, z, z_ind, z2

    end

    methods (Abstract)
        ind = get_ind(o);
        y = apply(o, data, params, c_vol)       
    end

    methods

        function o = mec_points(data)
            o.data = data;
            o.sz = data.dim(1:3);
            o.n_vol = data.n_vol;
        end

        % temporary placeholder
        function n = get.n(o)

            n = 0;
            for c_vol = 1:o.n_vol
                n = n + numel(o.get_ind(c_vol));
            end
            
        end

        function c = get_coords(o)
            c = o.coords;
        end
                
        function c = ind2coords(o, ind)
            
            % Use centered coordinates
            nx = o.sz(1);
            ny = o.sz(2);
            nz = o.sz(3);

            c.x = rem(ind-1, nx) + 1;
            c.y = rem(floor((ind-1)/nx), ny) + 1;
            c.z = floor((ind-1)/(nx*ny)) + 1;

            c.z_ind = c.z;

            % center linear terms
            c.x = c.x - nx / 2;
            c.y = c.y - ny / 2;
            c.z = c.z - nz / 2;

            % Add quadratic terms
            c.x2 = c.x.^2;
            c.y2 = c.y.^2;
            c.z2 = c.z.^2;
            c.xy = c.x .* c.y;
            c.xz = c.x .* c.z;
            c.yz = c.y .* c.z;

            % Add cubic terms
            c.x3   = c.x.^3;
            c.y3   = c.y.^3;
            c.z3   = c.z.^3;
            c.x2y  = c.x2 .* c.y;
            c.x2z  = c.x2 .* c.z;
            c.y2x  = c.y2 .* c.x;
            c.y2z  = c.y2 .* c.z;
            c.z2x  = c.z2 .* c.x;
            c.z2y  = c.z2 .* c.y;
            c.xyz  = c.x .* c.y .* c.z;
           
        end

        function c = distort(o, params, c_vol)

            c = o.get_coords();

            if (~isempty(params))
                c = params.distort(c, c_vol);
            end

            c.x = c.x + o.sz(1) / 2;
            c.y = c.y + o.sz(2) / 2;
            c.z = c.z + o.sz(3) / 2;            
            
        end

        function out = apply_one_vol(o, data, c_vol, params)

            if (nargin < 4), params = []; end % allows non-distortion

            c = o.distort(params, c_vol);
            
            out = data.interpolate(c.x, c.y, c.z, c_vol);

        end


    end

end