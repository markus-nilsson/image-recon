function out = imop_maxfilt3(data, windowSize)

assert(numel(windowSize) == 3, 'windowSize must be a 3-element vector');
fun = @(x) max(x(:));
out = nlfilter3(data, windowSize, fun);

end


function data_out = nlfilter3(data, win, fun)

% Custom 3D version of nlfilter
out = size(data);
pad = floor(win / 2);
data = padarray(data, pad, 'replicate');

data_out = zeros(out);
for c = 1:out(4)
    for x = 1:out(1)
        for y = 1:out(2)
            for z = 1:out(3)
                block = data(x:x+win(1)-1, y:y+win(2)-1, z:z+win(3)-1,c);
                data_out(x,y,z,c) = fun(block);
            end
        end
    end
end

end