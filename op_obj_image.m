classdef op_obj_image < op_obj

    properties (Access = protected)
        h_lhs;
        h_rhs;

        d_type = 'do_w_image_vector'; % applies to this data type

        aM; % adjusted operator
        aMT; % adjusted adjoint operator
        
    end

    methods

        function O = op_obj_image(S, h_rhs, h_lhs, n_k)

            if (nargin < 1), S = []; end
            if (nargin < 2), h_rhs = []; end
            if (nargin < 3), h_lhs = []; end
            if (nargin < 4), n_k = []; end

            O = O@op_obj();

            O.n_i = size(S,1);
            O.n_j = size(S,2);
            O.n_k = n_k;

            if (~isempty(h_lhs))
                O.n_l = h_lhs.dim(5);
            else
                O.n_l = [];
            end

            O.h_lhs = h_lhs;
            O.h_rhs = h_rhs;

            O.aM = S;
            O.aMT = S';

        end

        function data_rhs = init_x(O, a, b)
            if (nargin < 2), a = O.n_j; end
            if (nargin < 3), b = O.n_k; end
            h = O.h_rhs;
            h.dim(5) = b;
            data_rhs = do_w_image_vector(zeros(a,b), h);
        end

        function data_lhs = apply(O, data_rhs, ind)
            if (nargin < 3) || (isempty(ind)), ind = 1:size(data_rhs.w,2); end
            assert(my_isa(data_rhs, O.d_type), 'rhs not a %s', O.d_type);
            data_lhs = do_w_image_vector(O.i_apply(data_rhs.w(:,ind)), O.h_lhs);
        end

        function w = i_apply(O, w)
            w = O.aM * w;
        end

        function data_rhs = apply_adjoint(O, data_lhs, ind)
            if (nargin < 3) || (isempty(ind)), ind = 1:size(data_lhs.w,2); end
            assert(my_isa(data_lhs, O.d_type), 'rhs not a %s', O.d_type);
            data_rhs = do_w_image_vector(O.i_apply_adjoint(data_lhs.w(:,ind)), O.h_rhs);
        end

        function w = i_apply_adjoint(O, w)
            w = O.aMT * w;
        end

    end


    methods (Static)


        function ind = h2ind(h) % outputs coord in 4 x N
            
            dim = h.dim(2:4);

            x = repmat( permute((1:dim(1))', [1 2 3]), 1, dim(2), dim(3));
            y = repmat( permute((1:dim(2))', [2 1 3]), dim(1), 1, dim(3));
            z = repmat( permute((1:dim(3))', [2 3 1]), dim(1), dim(2), 1);

            ind = cat(2, x(:), y(:), z(:), ones(size(x(:))))';
        end

        function cm = h2cm(h) % header to camera matrix
            cm = [h.srow_x(:)'; h.srow_y(:)'; h.srow_z(:)'; 0 0 0 1];
        end

        function pc = ind2pc(ind, h) % index to physical coordinates
            C = op_obj_image.h2cm(h);
            ind(1:3,:) = ind(1:3,:) - 1;
            pc = C * ind;
        end

        function ind = pc2ind(pc, h) % physcial coordinates to indices
            C = op_obj_image.h2cm(h);
            ind = inv(C) * pc; %#ok<MINV> 
            ind(1:3,:) = ind(1:3,:) + 1;
        end

        function S = linear_interpolation(ind_lr, hr_dim, lr_dim)
            % ind_lr - 4 x N expected

            assert(size(ind_lr,1) == 4, 'size of ind_lr expected to be 4xN');

            ind_lr = ind_lr(1:3,:)';

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

            % shift matrix to go with the weight matrixs
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

    end
end
