classdef op_obj_image < op_obj

    properties (Access = protected)

        d_type = 'do_w_image_vector'; % applies to this data type

        aM; % adjusted operator
        aMT; % adjusted adjoint operator

    end

    properties

        S;

    end

    methods

        function O = op_obj_image(S, h_rhs, h_lhs, n_k)

            if (nargin < 1), S = []; end
            if (nargin < 2), h_rhs = []; end
            if (nargin < 3), h_lhs = []; end
            if (nargin < 4), n_k = []; end

            O = O@op_obj();

            if (~isempty(h_lhs)), O.n_i = prod(h_lhs.dim(2:4)); end
            if (~isempty(h_rhs)), O.n_j = prod(h_rhs.dim(2:4)); end
            O.n_k = n_k;

            if (~isempty(S))
                assert(O.n_i == size(S,1), 'wrong size');
                assert(O.n_j == size(S,2), 'wrong size');
            end

            if (~isempty(h_lhs))
                O.n_l = h_lhs.dim(5);
            else
                O.n_l = [];
            end

            O.h_lhs = h_lhs;
            O.h_rhs = h_rhs;

            if (canUseGPU)
                f = @(x) gpuArray(x);
            else
                f = @(x) x;
            end


            O.aM = f(S);
            O.aMT = f(S'); % wasteful on memory, but let it be for now

        end

        function data_lhs = apply(O, data_rhs, ind)

            if (nargin < 3), ind = []; end

            if (data_rhs.c_type == 2) % do_c; not pretty, redo (xxx) 

                data_lhs = do_c(numel(data_rhs));
                for c = 1:numel(data_rhs)
                    data_lhs.data_obj{c} = O.apply(data_rhs.data_obj{c}, ind);
                end

            elseif (data_rhs.c_type == 1) % do_w; not pretty, redo (xxx)

                if (isempty(ind)), ind = 1:size(data_rhs.w,2); end
                data_lhs = do_w_image_vector(O.i_apply(data_rhs, ind), O.h_lhs);

            else

                error('bad data type);')

            end

        end

        function w = i_apply(O, d, ind)
            w = O.aM * d.w(:, ind);
        end

        function data_rhs = apply_adjoint(O, data_lhs, ind)
            if (nargin < 3), ind = []; end
            
            if (data_lhs.c_type == 2) % do_c; not pretty, redo (xxx) 

                data_rhs = do_c(numel(data_lhs));
                for c = 1:numel(data_lhs)
                    data_rhs.data_obj{c} = O.apply_adjoint(data_lhs.data_obj{c}, ind);
                end
                
            
            elseif (data_lhs.c_type == 1) % do_w; not pretty, redo (xxx)

                if (isempty(ind)), ind = 1:size(data_lhs.w,2); end
            
                data_rhs = do_w_image_vector(O.i_apply_adjoint(data_lhs, ind), O.h_rhs);
            
            else
                error('bad data type');
            end

        end

        function w = i_apply_adjoint(O, d, ind)
            w = O.aMT * d.w(:, ind);
        end


        function S = get.S(O)
            S = O.aM;
        end


    end


    methods (Static)


        function ind = h2ic(h, ss) % outputs image coord in 4 x N

            if (nargin < 2), ss = 1; end
            
            dim = h.dim(2:4);

            x = permute((1:ss:dim(1))', [1 2 3]);
            y = permute((1:ss:dim(2))', [2 1 3]);
            z = permute((1:ss:dim(3))', [2 3 1]);

            dim = [numel(x) numel(y) numel(z)];

            x = repmat( x, 1, dim(2), dim(3));
            y = repmat( y, dim(1), 1, dim(3));
            z = repmat( z, dim(1), dim(2), 1);

            ind = cat(2, x(:), y(:), z(:), ones(size(x(:))));
        end

        function cm = h2cm(h) % header to camera matrix

            if (h.sform_code ~= 0)
                cm = op_obj_image.srow2cm(h);
            else
                cm = op_obj_image.quaternion2cm(h);
            end
        end

        function pc = ind2pc(ind, h) % index to physical coordinates
            C = op_obj_image.h2cm(h);
            ind(:,1:3) = ind(:,1:3) - 1;
            pc = ind * C';
        end

        function ind = pc2ind(pc, h) % physcial coordinates to indices
            C = op_obj_image.h2cm(h);
            ind = pc * inv(C)'; %#ok<MINV> 
            ind(:,1:3) = ind(:,1:3) + 1;
        end

        function S = linear_interpolation(ind_lr, hr_dim, lr_dim)
            % ind_lr - 4 x N expected

            assert(size(ind_lr,2) == 4, 'size of ind_lr expected to be 4xN');

            ind_lr = ind_lr(:,1:3);

            i = floor(ind_lr);
            w = ind_lr - i;

            % weight matrix for linear interpolation
            w = cat(2, ...
                (1-w(:,1)).*(1-w(:,2)).*(1-w(:,3)), ...
                ( w(:,1) ).*(1-w(:,2)).*(1-w(:,3)), ...
                (1-w(:,1)).*( w(:,2) ).*(1-w(:,3)), ...
                ( w(:,1) ).*( w(:,2) ).*(1-w(:,3)), ...
                (1-w(:,1)).*(1-w(:,2)).*( w(:,3) ), ...
                ( w(:,1) ).*(1-w(:,2)).*( w(:,3) ), ...
                (1-w(:,1)).*( w(:,2) ).*( w(:,3) ), ...
                ( w(:,1) ).*( w(:,2) ).*( w(:,3) ));

            % shift matrix to go with the weight matrix
            s = [...
                0 0 0;
                1 0 0;
                0 1 0;
                1 1 0;
                0 0 1;
                1 0 1;
                0 1 1;
                1 1 1]';

            n_lr = prod(lr_dim);
            n_hr = prod(hr_dim);

            S = spalloc(n_lr, n_hr, 8*n_hr);
            for c = 1:8

                i_lr_tmp = i + s(:,c)';

                ind_hr = ...
                    (i_lr_tmp(:,1) >= 1) & ...
                    (i_lr_tmp(:,2) >= 1) & ...
                    (i_lr_tmp(:,3) >= 1) & ...
                    (i_lr_tmp(:,1) <= lr_dim(1)) & ...
                    (i_lr_tmp(:,2) <= lr_dim(2)) & ...
                    (i_lr_tmp(:,3) <= lr_dim(3));

                ind_lr = ...
                    (i_lr_tmp(:,1)-1) + ...
                    (i_lr_tmp(:,2)-1) * lr_dim(1) + ...
                    (i_lr_tmp(:,3)-1) * lr_dim(2) * lr_dim(1) + ...
                    1;

                TMP = sparse(ind_lr(ind_hr), find(ind_hr), w(ind_hr,c), n_lr, n_hr);

                S = S + TMP;

            end

        end


        function R = quaternion2rotmat(h)
            % from https://github.com/bjw34032/oro.nifti/blob/master/R/quaternion.R

            assert(h.qform_code > 0, 'qform code must be above zero'); 

            tol = 1e-6;

            b = h.quatern_b;
            c = h.quatern_c;
            d = h.quatern_d;

            a = 1 - (b*b + c*c + d*d);

            if (a < tol)
                a = 1 / sqrt(b*b + c*c + d*d);
                b = a * b;
                c = a * c;
                d = a * d;
                a = 0;
            else
                a = sqrt(a);
            end

            R = [...
                a*a+b*b-c*c-d*d, 2*b*c+2*a*d, 2*b*d-2*a*c
                2*b*c-2*a*d, a*a+c*c-b*b-d*d, 2*c*d+2*a*b
                2*b*d+2*a*c, 2*c*d-2*a*b, a*a+d*d-c*c-b*b]';
        end

        function cm = quaternion2cm(h)

            qx = h.qoffset_x;
            qy = h.qoffset_y;
            qz = h.qoffset_z;

            dx = h.pixdim(2);
            dy = h.pixdim(3);
            dz = h.pixdim(4);


            qfac = h.pixdim(1);

            R = op_obj_image.quaternion2rotmat(h);

            assert(dx > 0, 'unexpected, negative pixel size');
            assert(dy > 0, 'unexpected, negative pixel size');
            assert(dz > 0, 'unexpected, negative pixel size');

            c_case = 1;
            switch (c_case) 
                case 1
                    if (qfac < 0), dz = -dz; end
                    R = R * diag([dx dy dz]);
                case 2
                    if (qfac < 0), dz = -dz; end
                    R(:,1) = R(:,1) * dx;
                    R(:,2) = R(:,2) * dy;
                    R(:,3) = R(:,3) * dz;
                    
            end

            cm = [R [qx qy qz]'; 0 0 0 1];

        end

        function cm = srow2cm(h)
            cm = [h.srow_x(:)'; h.srow_y(:)'; h.srow_z(:)'; 0 0 0 1];
        end

    end
end
