classdef mec_params_global_orders < mec_params_base

    properties
        ec_order = 0;
    end

    methods

        function o = mec_params_global_orders(data, ec_order)
            o = o@mec_params_base(data);
            o.ec_order = ec_order;

            o = o.reset();

            o.opt_opts.Display = 'iter';
            o.opt_opts.TypicalX = ones(size(o.get_x(1)));
            
        end

        function o = reset(o)
        
            switch (o.ec_order)
                case 0
                    o = o.update(zeros(1,3));
                case 1
                    o = o.update(zeros(1,12));
                case 2
                    o = o.update(zeros(1,30));
                case 3
                    o = o.update(zeros(1,60));
                otherwise
                    error('undefined');
            end
            o.enabled = false;
        end        

        function o = update(o, x, ~)
            o.x = x;
            o.enabled = true;
        end

        function x = get_x(o, ~)
            x = o.x;
        end

        function p = get_p(o, ~)

            x = o.get_x(o);

            % improve scaling, try to keep normal values at 1
            p = x;

            if (o.ec_order >= 0) % 0th order
                p(1:3) = p(1:3) * 1e0; 
            end

            if (o.ec_order >= 1) % 1st order
                p(4:12) = p(4:12) * 1e-2; % 1st order
            end

            if (o.ec_order >= 2)
                p(13:30) = p(13:30) * 1e-3; % 2nd order terms
            end

            if (o.ec_order >= 3)
                p(31:60) = p(31:60) * 1e-4; % 2nd order terms
            end
            
        end

        function [d,i] = distort(o, c, c_vol)

            if (~isempty(o.prior_param))
                c = o.prior_param.distort(c, c_vol);
            end

            d = c;

            if (~o.enabled), return; end
           
            p = o.get_p();
            u = o.u;
            i = 0;

            % Homogeneous eddy currents
            d0 = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;

            d.y = d.y + d0; % homogeneous eddies

            if (o.ec_order == 0), return; end

            % Linear eddy currents
            gx = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            gy = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            gz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
                        
            d.y = d.y + ...
                gx .* c.x + ... % shear
                gy .* c.y + ... % scale
                gz .* c.z;      % tilt

            if (o.ec_order == 1), return; end
            

            % Second order terms
            g_xx = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_xy = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_xz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_yy = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_yz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_zz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;

            d.y = d.y ...
                + g_xx .* c.x2 ...
                + g_yy .* c.y2 ...
                + g_zz .* c.z2 ...
                + g_xy .* c.xy ...
                + g_xz .* c.xz ...
                + g_yz .* c.yz;

            if (o.ec_order == 2), return; end    


            % Third order terms
            g_xxx = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_yyy = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_zzz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;

            g_xxy = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_xxz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_yyx = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_yyz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_zzx = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;
            g_zzy = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;

            g_xyz = u(c_vol, :) * p( (1:3) + i)'; i = i + 3;

            d.y = d.y ...
                + g_xxx .* c.x3 ...
                + g_yyy .* c.y3 ...
                + g_zzz .* c.z3 ...
                + g_xxy .* c.x2y ...
                + g_xxz .* c.x2z ...
                + g_yyx .* c.y2x ...
                + g_yyz .* c.y2z ...
                + g_zzx .* c.z2x ...
                + g_zzy .* c.z2y ...
                + g_xyz .* c.xyz;

            if (o.ec_order == 3), return; end            
            
        end

        function o = optimize(o, data_mov, data_ref, points)

            % 1. Grab reference values 
            y = points.apply(data_ref);

            % Define objective function
            f_obj = @(x) o.opt_obj(x, y, data_mov, points, []);

            % Rescale it
            sc = o.opt_rescale(f_obj, o.get_x());
            f_obj = @(x) f_obj(x) / sc;

            % 2. Optimize
            tic;
            x_hat = lsqnonlin(f_obj, o.get_x(), [], [], o.opt_opts);
            t = toc;

            % 3. Display
            o.opt_disp(f_obj, x_hat, t);

            % 4. Store
            o = o.update(x_hat);

        end        

        function plot(o)

            msf_clf;

            p_tmp = o.get_p();
            p = zeros(1, 60);
            p(1:numel(p_tmp)) = p_tmp;

            % 
            plot_shift = [1  8 9 10  15:20];

            yscl = [1   0.02 0.02 0.02  1e-3 1e-3 1e-3  1e-3 1e-3 1e-3 ];

            t_str = {'dy', ...
                'shear', 'scale', 'trans', ...
                'xx', 'xy', 'xz', 'xy', 'xz', 'yz'};

            for c = 1:10

                subplot(3,7, plot_shift(c));

                plot(1:3, [0 0 0], 'k:', 1:3, p( (1:3) + (c-1) * 3), 'o-');
                ylim(yscl(c) * [-1 1]);
                xlim([0.5 3.5]);
                set(gca,'xtick', 1:3, 'xticklabel', {'x','y','z'});
                title({t_str{c}, ''});

            end

            pause(0.05);

        end

    end


    methods (Static)

        function stop = my_callback(p, x, optimValues, state)

            stop = 0;

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