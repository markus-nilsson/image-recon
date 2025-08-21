classdef mec_params_global_orders < mec_params_base
    %
    % MEC (motion/eddy current) global-parameter distortion model.
    %
    % Purpose
    %   Encapsulates a global polynomial model for through-phase-encode
    %   displacement d.y induced by eddy currents, parameterized per diffusion
    %   encoding (row in o.u). Supports 0th–3rd order spatial polynomials in
    %   x,y,z. Parameters are estimated globally via the base-class optimizer.
    %
    % Coordinate conventions / inputs
    %   - c: struct of sample coordinates and precomputed monomials
    %       required fields (always):
    %         c.x, c.y, c.z, c.z_ind (slice index),
    %       for ≥2nd order: c.x2,c.y2,c.z2,c.xy,c.xz,c.yz,
    %       for ≥3rd order: c.x3,c.y3,c.z3,c.x2y,c.x2z,c.y2x,c.y2z,c.z2x,c.z2y,c.xyz.
    %     Coordinates should be centered (zero-mean); see calling code.
    %   - u: (nvol × 3) diffusion-encoding table held in o.u; each row is the
    %     (x,y,z) component of the unit gradient (or any 3-vector “encoding”).
    %     All parameter blocks are modulated by u for the current volume.
    %
    % Parameterization
    %   Order 0 (3 params):      [d0x d0y d0z] — homogeneous shift per encoding
    %   Order 1 (9 params):      [gx gy gz] for coefficients multiplying x,y,z
    %                            giving shear/scale/tilt along phase-encode (y)
    %   Order 2 (18 params):     [g_xx g_xy g_xz g_yy g_yz g_zz]
    %   Order 3 (30 params):     [g_xxx g_yyy g_zzz g_xxy g_xxz g_yyx g_yyz g_zzx g_zzy g_xyz]
    %   Total per order: 0→3, 1→12, 2→30, 3→60 parameters.
    %
    % Scaling
    %   get_p() rescales raw x to p to keep typical magnitudes ≈1 across
    %   orders (eases optimization). Tune the 1e0/1e-2/1e-3/1e-4 as needed.
    %
    % Output of distort()
    %   Returns d struct with updated coordinates (only d.y is modified) and an
    %   index i used internally to walk the parameter vector.
    %
    % Notes
    %   - This class only models global, encoding-modulated displacements in
    %     the PE direction. Fieldmaps or local models live elsewhere.
    %   - The optimizer callback animates parameter evolution.
    %   - Keep ec_order consistent with the size of x in reset().
    %
    properties
        ec_order = 0;   % maximum polynomial order (0..3)
    end

    methods

        function o = mec_params_global_orders(data, ec_order)
            % Constructor: set base, order, default state, and optimizer opts
            o = o@mec_params_base(data);
            o.ec_order = ec_order;

            o = o.reset();             % allocate x with correct length

            o.opt_global = 1;          % signal global optimization block

            % Optimizer housekeeping
            o.opt_opts.Display  = 'iter';
            o.opt_opts.TypicalX = ones(size(o.get_x(1))); % scale hint
        end

        function o = reset(o)
            % Allocate parameter vector x according to ec_order, then disable
            % the block until explicitly updated/used.
            switch (o.ec_order)
                case 0, o = o.update(zeros(1,  3));
                case 1, o = o.update(zeros(1, 12));
                case 2, o = o.update(zeros(1, 30));
                case 3, o = o.update(zeros(1, 60));
                otherwise, error('undefined');
            end
            o.enabled = false;
        end

        function o = update(o, x, ~)
            % Set raw parameter vector and enable the block
            o.x = x;
            o.enabled = true;
        end

        function x = get_x(o, ~)
            % Return raw (unscaled) parameter vector
            x = o.x;
        end

        function p = get_p(o, ~)
            % Return scaled parameter vector used internally in distort()
            % Scaling aims for order-balanced magnitudes to improve conditioning
            x = o.get_x(o);
            p = x; % start from raw

            if (o.ec_order >= 0) % 0th order (homogeneous)
                p(1:3) = p(1:3) * 1e0;  % typical ~1
            end
            if (o.ec_order >= 1) % 1st order (linear terms)
                p(4:12) = p(4:12) * 1e-2;
            end
            if (o.ec_order >= 2) % 2nd order (quadratic terms)
                p(13:30) = p(13:30) * 1e-3;
            end
            if (o.ec_order >= 3) % 3rd order (cubic terms)
                p(31:60) = p(31:60) * 1e-4;
            end
        end

        function [d,i] = distort(o, c, c_vol)
            % Compute PE-direction displacement for volume c_vol
            % Inputs
            %   c: coordinate/monomial struct (see header)
            %   c_vol: index into o.u selecting current encoding row
            % Outputs
            %   d: struct with updated coordinates (d.y modified in-place)
            %   i: last used index into p (mostly for debugging/inspection)

            % Pass-through of positions (only y is distorted here)
            d.x = c.x; d.y = c.y; d.z = c.z; d.z_ind = c.z_ind;
            if (~o.enabled), return; end

            p = o.get_p();  % scaled param vector
            u = o.u;        % encoding table (nvol×3)
            i = 0;          % running index into p

            % -----------------------------
            % 0th order: homogeneous shift (per encoding)
            % d0 = dot(u, p0)
            d0 = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            d.y = d.y + d0;
            if (o.ec_order == 0), return; end

            % -----------------------------
            % 1st order: linear terms — shear/scale/tilt
            % gx multiplies x, gy multiplies y, gz multiplies z
            gx = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            gy = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            gz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            d.y = d.y + gx.*c.x + gy.*c.y + gz.*c.z;
            if (o.ec_order == 1), return; end

            % -----------------------------
            % 2nd order: quadratic terms
            g_xx = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_xy = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_xz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_yy = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_yz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_zz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            d.y = d.y + g_xx.*c.x2 + g_yy.*c.y2 + g_zz.*c.z2 ...
                + g_xy.*c.xy + g_xz.*c.xz + g_yz.*c.yz;
            if (o.ec_order == 2), return; end

            % -----------------------------
            % 3rd order: cubic terms
            g_xxx = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_yyy = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_zzz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_xxy = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_xxz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_yyx = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_yyz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_zzx = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_zzy = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            g_xyz = u(c_vol, :) * p((1:3)+i)'; i = i + 3;
            d.y = d.y + g_xxx.*c.x3 + g_yyy.*c.y3 + g_zzz.*c.z3 ...
                + g_xxy.*c.x2y + g_xxz.*c.x2z ...
                + g_yyx.*c.y2x + g_yyz.*c.y2z ...
                + g_zzx.*c.z2x + g_zzy.*c.z2y ...
                + g_xyz.*c.xyz;
            if (o.ec_order == 3), return; end
        end

        function plot(o)
            % Quick-look diagnostic plot of parameter blocks vs encoding axes
            msf_clf;

            p_tmp = o.get_p();
            p = zeros(1,60); p(1:numel(p_tmp)) = p_tmp; % pad for plotting

            % panel positions; hand-tuned indices
            plot_shift = [1 8 9 10 15:20];

            % y-limits per panel (rough scales by order)
            yscl = [1 0.04 0.04 0.04 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3];

            t_str = {'dy', 'shear','scale','trans', ...
                'xx','xy','xz','yy','yz','zz'};

            for c = 1:10
                subplot(3,7, plot_shift(c));
                plot(1:3, [0 0 0], 'k:', 1:3, p((1:3)+(c-1)*3), 'o-');
                ylim(yscl(c) * [-1 1]); xlim([0.5 3.5]);
                set(gca,'xtick',1:3,'xticklabel',{'x','y','z'});
                title({t_str{c}, ''});
            end
            pause(0.05);
        end
    end

    methods (Static)
        function stop = my_callback(p, x, optimValues, state)
            % Optimizer callback: update live plot every iteration
            stop = 0; %#ok<NASGU>
            switch (state)
                case 'init'
                    return;
            end
            figure(my_fig('mec_global'));
            p.update(x).plot();
            pause(0.05);
        end
    end
end
