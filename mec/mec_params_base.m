classdef mec_params_base

    properties
        x;    % parameters to fit
        data; % handle to data
        
        % Current and previous gradient directions
        u;
        v;
        w;

        enabled = false 

        % Optimization related 
        opt_opts = optimoptions('lsqnonlin', ...
            'Algorithm','levenberg-marquardt', ...
            'Display','off', ...
            'FunctionTolerance',1e-5, ...
            'StepTolerance',1e-5, ...
            'MaxIterations',40);

        opt_sc = 1; % optimization scaling

        opt_global = false;
    end
    
    methods

        function o = mec_params_base(data)

            % Store data
            o.data = data;

            if (~isfield(o.data.xps, 'u_from_bvec'))
                warning('check your xps, u_from_bvec missing');
            else

                if (sum( (o.data.xps.u(:) - o.data.xps.u_from_bvec(:)).^2 ) > 0.001)
                    error('xps.u corrupted, fixing it')
                end

            end

            % Prepare gradient direction sets
            o.u = o.data.xps.u .* ...
                repmat(sqrt(o.data.xps.b / max(o.data.xps.b)), 1, 3);

            o.v = circshift(o.u, 1, 1); o.v(1,:) = 0;
            o.w = circshift(o.v, 1, 1); o.w(1,:) = 0;
            
            % Reset parameter                        
            o = o.reset();

            % Tell optimizer about typical x
            o.opt_opts.TypicalX = ones(size(o.get_x(1)));
            
        end

    end

    methods (Abstract)
        reset(o) % null parameters
        update(o, x, c_vol) % store fit result
        get_x(o, c_vol); % fit parameters
        get_p(o, c_vol) % rescale from fit variables
        distort(o, x, y, z, c_vol) % distort volume
        plot(o) 
    end

    methods % for optimization

        function o = optimize(o, data_mov, data_ref, points)

            if (o.opt_global)
                c_vols = 0;
            else
                c_vols = 1:data_mov.n_vol;
            end

            for c_vol = c_vols

                % 1. Grab reference values and rescale obj function
                y = points.apply(data_ref, [], c_vol);

                % Define objective function
                f_obj = @(x) o.opt_obj(x, y, data_mov, points, c_vol);

                % Rescale it
                sc = o.opt_rescale(f_obj, o.get_x(c_vol));
                f_obj = @(x) f_obj(x) / sc;

                % 2. Optimize
                tic;
                x_hat = lsqnonlin(f_obj, o.get_x(c_vol), [], [], o.opt_opts);
                t = toc;

                % 3. Display
                o.opt_disp(f_obj, x_hat, t, c_vol);

                % 4. Store
                o = o.update(x_hat, c_vol);                

            end
        end

        function L = opt_obj(o, x, y, data_mov, points, c_vol)

            y_hat = points.apply(data_mov, o.update(x, c_vol), c_vol);

            L = y - y_hat;

        end

        function sc = opt_rescale(~, f_obj, x0)
            sc = norm(f_obj(x0));
            sc = max(0.001, sc);
        end

        function opt_disp(o, f_obj, x_hat, t, c_vol)

            if (nargin < 5), c_vol = []; end

            x0 = o.get_x(c_vol);

            disp(sprintf('%i: %1.5f, %1.2f in %1.1f s', ...
                c_vol, ...
                norm(f_obj(x0)).^2, ...
                (norm(f_obj(x_hat)) / norm(f_obj(x0))).^2, ...
                t));

        end

    end

end
