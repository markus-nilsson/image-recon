classdef data_obj_image_cell < data_obj_cell

    properties
        h
    end

    methods

        function h = get.h(O)
            h = cell(1, numel(O.data_obj));
            for c = 1:numel(h)
                h{c} = O.data_obj{c}.h;
            end
        end

    end

end