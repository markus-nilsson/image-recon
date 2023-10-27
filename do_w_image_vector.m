classdef do_w_image_vector < do_w

    properties (GetAccess = public, SetAccess = protected)
        h;
    end

    methods

        function O = do_w_image_vector(w, h)
            assert(size(w,1) == prod(h.dim(2:4)), 'wrong size');
            O = O@do_w(w);            
            h.dim(5) = size(w,2);
            O.h = h;            
        end

        function O_new = new(O, w)
            assert(prod(O.h.dim(2:4)) == size(w, 1), 'size error');
            O.h.dim(5) = size(w,2);
            O_new = do_w_image_vector(w, O.h);
        end

        function O = add(O, w, l)
            if (my_isa(w, 'do_w')), f = @(x) x.w; else, f = @(x) x; end
            O.w(:,l) = O.w(:,l) + f(w);
        end

        function I = imreshape(O)
            I = reshape(O.w, [O.h.dim(2:4)' size(O.w,2)]);            
        end

        function O_new = zeros(O)
            O_new = do_w_image_vector(zeros(size(O.w)), O.h);
        end     

        function d = dim(O, ind)
            d = O.h.dim(2:5)';
            if (nargin > 1)
                d = d(ind);
            end
        end

        function d = pixdim(O, ind)
            d = O.h.pixdim(2:5)';
            if (nargin > 1)
                d = d(ind);
            end
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
        function do = new_hr_from_lr(data_lr)

            h_hr = data_lr.h;
            tmp = single(...
                h_hr.dim(4))*h_hr.pixdim(4)/(h_hr.pixdim(2));
            h_hr.dim(4) = int16(tmp);
            h_hr.pixdim(4) = data_lr.h.pixdim(2);

            I = zeros(double(h_hr.dim(2:5)'));

            do = do_w_image_volume(I,h_hr);

        end
    end
end

