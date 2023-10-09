classdef do_w_image_volume < do_w_image_vector
    
    methods

        function O = do_w_image_volume(I, h)

            I = reshape(I, [prod(size(I,1,2,3)) size(I,4)]);
            O = O@do_w_image_vector(I, h);

        end

        function O = trim(O, ir, jr, kr, vr)

            if (nargin < 2) || (isempty(ir)), ir = 1:O.h.dim(2); end
            if (nargin < 3) || (isempty(jr)), jr = 1:O.h.dim(3); end
            if (nargin < 4) || (isempty(kr)), kr = 1:O.h.dim(4); end
            if (nargin < 5) || (isempty(vr)), vr = 1:O.h.dim(5); end

            I = O.imreshape();
            I = I(ir, jr, kr, vr);
            O.w = reshape(I, [prod(size(I,1,2,3)) size(I,4)]);
            O.n_vox = size(O.w, 1);

            O.h.dim(2:5) = size(I, 1, 2, 3, 4);

            % adjust offsets
            ind = [min(ir)-1 min(jr)-1 min(kr)-1 1]';
            cm = [O.h.srow_x'; O.h.srow_y'; O.h.srow_z';  0 0 0 1];
            physc = cm * ind;

            O.h.qform_code = 0;
            O.h.srow_x(4) = physc(1);
            O.h.srow_y(4) = physc(2);
            O.h.srow_z(4) = physc(3);
            
        end

    end

    methods (Static)

        % not ideal place for this
        function data_obj = new_hr_from_lr(data_lr)

            h_hr = data_lr.h;
            tmp = single(...
                h_hr.dim(4))*h_hr.pixdim(4)/(h_hr.pixdim(2));
            h_hr.dim(4) = int16(tmp);
            h_hr.pixdim(4) = data_lr.h.pixdim(2);

            % also adjust srow's
            h_hr.srow_z(3) = h_hr.pixdim(4);

            X = [h_hr.srow_x'; h_hr.srow_y'; h_hr.srow_z';];
            X = X(1:3, 1:3);
            X = X - diag(diag(X));

            if (sum(X(:).^2) > 1e-5)
                error('need transverse slices for now');
            end


            I = zeros(double(h_hr.dim(2:5)'));

            data_obj = do_w_image_volume(I,h_hr);

        end
    end
end

