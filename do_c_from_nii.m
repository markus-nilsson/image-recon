classdef do_c_from_nii < do_c

    properties
        fns = {};
        h;
    end

    methods

        function O = do_c_from_nii(fns)
            assert(iscell(fns), 'input must be cells with file names');

            O = O@do_c(numel(fns));
            O.fns = fns;

            for c = 1:numel(fns)
                O.data_obj{c} = do_w_from_nii(fns{c});
            end

        end

        function O = trim(O, ir, jr, kr, vr)
            for c = 1:numel(O)
                O.data_obj{c}.trim(ir, jr, kr, vr);
            end
        end

        function h = get.h(O)
            h = cell(1, numel(O));
            for c = 1:numel(h)
                h{c} = O.data_obj{c}.h;
            end
        end

    end
end