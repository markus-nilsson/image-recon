classdef op_obj_smr < op_obj_dict

    properties % this need to be fixed up
      
        M_odf;
        DBmerge;

        k02_ind;
        odf_ind;

        c_type;

        mp;

        CEk_full;

    end
    
    methods
        
        function O = op_obj_smr(xps, n_vox)
            % function O = op_obj_smr(xps, n_vox)

            if (nargin < 2), n_vox = 1; end

            % compute basis functions for stick and zeppelin separately
            
            n = 8; % increase in later experiments

            f_list = linspace(0, 1, n);
            di_s_list = linspace(0.2e-9, 1.3e-9, n);
            di_z_list = linspace(0.0e-9, 3.3e-9, n);
            dd_z_list = linspace(0.0, 0.8700, n);
            t2_s_list = linspace(30.e-3, 300e-3, n);
            t2_z_list = 1./linspace(1./30.e-3, 1./500e-3, n);

            disp('using focused set');
            di_s_list = linspace(0.6e-9, 1.0e-9, n);
            di_z_list = linspace(0.5e-9, 1.5e-9, n);
            dd_z_list = linspace(0.0, 0.8700, n);
            t2_s_list = linspace(30.e-3, 150e-3, n);
            t2_z_list = 1./linspace(1./60.e-3, 1./300e-3, n);
            

            % stick signal database
            S_s = zeros(2*xps.n, n, n);
            for i = 1:n
                m = [1 1 di_s_list(i) 3e-9 0 0 0 0 0 0 inf 80e-3];
                [~, export] = dtd_smr_1d_fit2data(m, xps, 1);

                tmp = cat(1, ...
                    exp(export.a_s) .* export.i0_s / sqrt(4*pi), ...
                    exp(export.a_s) .* export.i2_s);

                for j = 1:n
                    t2w = exp(-xps.te / t2_s_list(j));
                    S_s(:,i,j) = tmp .* cat(1, t2w, t2w);
                end
            end

            % zeppelin signal database
            S_z = zeros(2*xps.n, n, n, n);
            for i = 1:n
                for j = 1:n
                    m = [1 3 3e-9 di_z_list(i) dd_z_list(j) 0 0 0 0 0 80e-3 inf];
                    [~,export] = dtd_smr_1d_fit2data(m, xps, 1);

                    tmp = cat(1, ...
                        exp(export.a_z) .* export.i0_z / sqrt(4*pi), ...
                        exp(export.a_z) .* export.i2_z);

                    for k = 1:n
                        t2w = exp(-xps.te / t2_z_list(k));
                        S_z(:,i,j,k) = tmp .* cat(1, t2w, t2w); 
                            
                    end
                end
            end

            P_z = repmat(S_z,  [1 1 1 1 n n]);
            P_z = permute(P_z, [1 5 6 2 3 4]);

            P_s = repmat(S_s,  [1 1 1 n n n]);

            S_k = zeros(2*xps.n, n, n, n, n, n, n);
            for i = 1:numel(f_list)
                S_k(:, :, :, :, :, :, i) = ...
                    f_list(i) * P_s + (1-f_list(i)) * P_z;

            end

            
            mp_list = {...
                di_s_list, ...
                t2_s_list, ...
                di_z_list, ...
                dd_z_list, ...
                t2_z_list, ...
                f_list};
            





            % build the sph basis function
            
            % Convert gradient vectors from Cartesian to Spherical
            x       = xps.u(:,1);
            y       = xps.u(:,2);
            z       = xps.u(:,3);
            [phi, theta] = cart2sph(x, y, z);
            phi     = phi + pi;       % [-pi pi]  --> [0 2pi], longitude
            theta   = pi / 2 - theta; % latitude  --> co-latitude

            %
            y20 = dtd_smr_spha(2,  0, theta, phi);
            y21 = dtd_smr_spha(2,  1, theta, phi);
            y22 = dtd_smr_spha(2,  2, theta, phi);

            T = sqrt(4*pi/5) * cat(2, ...
                y20, ...
                ( 2) * real(y21), ...
                (-2) * imag(y21), ...
                ( 2) * real(y22), ...
                (-2) * imag(y22));

            T = cat(2, ...
                cat(1, ones(xps.n,1), zeros(xps.n,1) ), ...
                cat(1, zeros(size(T)), T) );


            O = O@op_obj_dict(S_k, mp_list, 35, n_vox);            
            
            O.M_odf = T;


            % control
            i1 = 4;
            i2 = 4;
            i3 = 4;
            i4 = 4;
            i5 = 4;
            i6 = 4;
            odf_pars = [0 0 0 0 0] + 0.5;
            mx = [1 f_list(i1) di_s_list(i2) di_z_list(i3) dd_z_list(i4) odf_pars t2_s_list(i5) t2_z_list(i6)];
            [s1,e] = dtd_smr_1d_fit2data(mx, xps, 1);
            s2 = S_k(:, i2, i5, i3, i4, i6, i1);


            % Build splitter matrix
            [~,c_list, id_ind] = mdm_pa_ind_from_xps(xps);

            % multiplication with Ma performs powder averaging, Mb restores
            Ma = zeros(numel(c_list), xps.n);
            Mb = zeros(numel(c_list), xps.n);

            for c = c_list'
                Ma(c, :) = (id_ind == c) / (sum(id_ind == c));
                Mb(c, :) = (id_ind == c);
            end

            M2 = Mb' * Ma;

            M3 = eye(size(Ma,2)) - M2;

            % decomposes the data into to parts: with and without null component
            % the transpose?
            M4 = cat(1, M2, M3);

            O.DBmerge = M4';

            % final logistics
            O.k02_ind = 1:size(O.M,2);
            O.odf_ind = size(O.M,2) + (1:6);

            O.c_type = 1;

            O.n_vox = n_vox;

            O.mp = O.reset_mp();
        end

        function mp = reset_mp(O)
            n_mp = size(O.M,2) + numel(O.odf_ind);
            mp = zeros(n_mp, O.n_vox);
            mp(O.odf_ind(1),:) = 1; % odf powder
        end


        function x = init_x(O, c_mode)
%             c_mode = 2;

            if (nargin < 2), c_mode = 1; end

            switch (c_mode)
                case 1
                    x = O.init_zero();
                case 2
                    x = repmat(mean(O.CEk,1)', [1 O.n_vox]);
                    x  = x .* mean(data ./ (O * x),1);                    
            end
        end

        function x = init_zero(O)
            switch (O.c_type)
                case 1
                    x = zeros(size(O.M,2), O.n_vox);
                case 2
                    x = zeros(size(O.M_odf,2), O.n_vox);
            end
        end      

        function y_dwi = apply(O,y_mp)

            switch (O.c_type)
                case 0
                    s_dwi   = O.M     * y_mp(O.k02_ind,:);
                    tmpodf  = O.M_odf * y_mp(O.odf_ind,:);

                case 1 % kernel
                    s_dwi   = O.M   * y_mp;
                    tmpodf  = O.M_odf * O.mp(O.odf_ind,:);

                case 2 % odf
                    s_dwi   = O.M     * O.mp(O.k02_ind,:);
                    tmpodf  = O.M_odf * y_mp;

            end

            y_dwi =  O.DBmerge * (s_dwi .* tmpodf);
        end
        

        function y_mp = apply_adjoint(O, y_dwi)

            switch (O.c_type)

                case 1
                    tmpodf = O.M_odf * O.mp(O.odf_ind,:);
                    y_mp   = O.M'  * (tmpodf .* (O.DBmerge' * y_dwi));

                case 2
                    s_dwi = O.M    * O.mp(O.k02_ind,:);
                    y_mp  = O.M_odf' * (s_dwi  .* (O.DBmerge' * y_dwi));

            end

        end

    end
end
