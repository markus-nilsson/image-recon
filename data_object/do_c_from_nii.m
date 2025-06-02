classdef do_c_from_nii < do_c

    properties
        fns = {};
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

        function O = powder_average(O)

            for c = 1:numel(O.data_obj)
                O.data_obj{c}.powder_average();
            end

        end

    end
end