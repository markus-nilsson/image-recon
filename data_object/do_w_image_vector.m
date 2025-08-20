classdef do_w_image_vector < do_w 

    properties (GetAccess = public, SetAccess = protected)
        h; % nifti header       
    end

    methods

        function O = do_w_image_vector(w, h, xps)
            if (nargin < 3), xps = []; end
            assert(size(w,1) == prod(h.dim(2:4)), 'wrong size');
            O = O@do_w(w, xps);            
            h.dim(5) = size(w,2);
            O.h = h;            
        end

        function O_new = new(O, w)
            if (ismatrix(w))
                assert(prod(O.h.dim(2:4)) == size(w, 1), 'size error');
            else
                assert(all(size(w, [1 2 3 4])' == O.h.dim(2:5)), 'size error');
                w = reshape(w, prod(size(w, [1 2 3])), size(w,4));
            end
            h = O.h;
            h.dim(5) = size(w,2);
            O_new = do_w_image_vector(w, h);
            O_new.xps = O.xps; % fix this
        end

        function O = add(O, w, l)
            if (isnumeric(w)), f = @(x) x; else f = @(x) x.w; end
            O.w(:,l) = O.w(:,l) + f(w);
        end

        function I = imreshape(O, c_vol)
            arguments
                O
                c_vol = 1:size(O.w,2)
            end

            I = gather(reshape(O.w(:,c_vol), [O.h.dim(2:4)' numel(c_vol)]));            
        end

        function O_new = zeros(O)
            O_new = do_w_image_vector(zeros(size(O.w)), O.h);
        end   

        function O = subsample(O, ind)

            O = subsample@do_w(O, ind);
            O.h.dim(5) = size(O.w, 2);
            
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

        function O = trim_copy(O, ir, jr, kr, vr)
            if (nargin < 5), vr = []; end
            O = trim(O, ir, jr, kr, vr, @(w) do_w_image_vector(w, O.h));
        end
        

        function O = trim(O, ir, jr, kr, vr, f)

            % keep track of input info, to accelerate trimming
            ii = 'yyyy'; 

            if (nargin < 2) || (isempty(ir)), ir = 1:O.h.dim(2); ii(1) = 'n'; end
            if (nargin < 3) || (isempty(jr)), jr = 1:O.h.dim(3); ii(2) = 'n'; end
            if (nargin < 4) || (isempty(kr)), kr = 1:O.h.dim(4); ii(3) = 'n'; end
            if (nargin < 5) || (isempty(vr)), vr = 1:O.h.dim(5); ii(4) = 'n'; end
            if (nargin < 6), f = @(w) O.set_w(w); end

            xps = O.xps;

            % Trim efficiently
            switch (ii) % implement [0 0 1 1] later
                case 'nnny'
                    O = f(O.w(:,vr));

                otherwise
                    I = O.imreshape();
                    I = I(ir, jr, kr, vr);
                    O = f(reshape(I, [prod(size(I,1,2,3)) size(I,4)]));

            end

            % Adjust info and header including offsets
            ind = [min(ir)-1 min(jr)-1 min(kr)-1 1]';
            cm = [O.h.srow_x'; O.h.srow_y'; O.h.srow_z';  0 0 0 1];
            physc = cm * ind;

            O.h.qform_code = 0;
            O.h.srow_x(4) = physc(1);
            O.h.srow_y(4) = physc(2);
            O.h.srow_z(4) = physc(3);
            O.h.dim(2:5) = [numel(ir) numel(jr) numel(kr) numel(vr)];
            

            if (~isempty(xps))

                if (numel(vr) ~= xps.n)
                    vr_ind = (1:xps.n) == 0;
                    for c = 1:numel(vr)
                        vr_ind(vr(c)) = 1 > 0;
                    end
                    
                    if (numel(vr_ind) ~= xps.n)
                        error('something badly specified');
                    end

                else
                    vr_ind = vr;
                end

                O.xps = mdm_xps_subsample(xps, vr_ind);
            end
            
        end       


        function w = patch(O, i,j,k,l,s)

            if (s > 0)
                ir = i + (-s:s)';
                jr = j + (-s:s)';
                kr = k + (-s:s)';
            else
                ir = i;
                jr = j;
                kr = k;
            end

            ir = ir( (ir >= 1) & (ir <= O.dim(1)));
            jr = jr( (jr >= 1) & (jr <= O.dim(2)));
            kr = kr( (kr >= 1) & (kr <= O.dim(3)));

            sub = combvec(ir', jr', kr');
            ind = sub2ind(O.dim(1:3), sub(1, :), sub(2,:), sub(3, :));

            w = reshape(O.w(ind, l), numel(ir), numel(jr), numel(kr));

        end

        function O = recast(O, target)

            O.w = cast(O.w, 'like', target.w);
            O.h.datatype = target.h.datatype;
            
        end

        function O = rescale(O, target)

            O.w = O.w / target.h.scl_slope * O.h.scl_slope;
            O.h.scl_slope = target.h.scl_slope;

        end


        function nii_fn = nii_write(obj, nii_fn, type)

            if (nargin < 3), type = []; end % cast
            if (isempty(type)), type = 'undefined'; end

            this_h = obj.h;
            this_I = obj.imreshape();

            switch (type)

                case 'single'
                    this_h.datatype = 16;
                    this_I = single(this_I);

                case 'undefined'
                    % do nothing

                otherwise
                    error('not implemented');
            end

            mdm_nii_write(this_I, nii_fn, this_h);

        end

        function show(O, c_dim, k, c_vol)

            arguments
                O
                c_dim = []
                k = []
                c_vol = []
            end

            if (isempty(c_vol)), c_vol = 1; end

            msf_imagesc(O.imreshape(c_vol), c_dim, k);
            % caxis([0 1600]);

        end

        function play(O, c_dim, k, n_rep, cscl)
        % function play(O, c_dim, k, n_rep, cscl)

            arguments
                O
                c_dim = []
                k = []
                n_rep = 2
                cscl = [];
            end

            for c_rep = 1:n_rep
                for c_vol = 1:size(O.w,2)
                    O.show(c_dim, k, c_vol);
                    title(num2str(c_vol));
                    if (~isempty(cscl)), caxis(cscl); end
                    pause(0.1);
                end
            end

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

