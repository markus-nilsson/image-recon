function y = imop_mean_closest_two(x, dim)


if nargin < 2
    dim = find(size(x) ~= 1, 1); % default to first non-singleton dimension
    if isempty(dim), dim = 1; end
end

if size(x, dim) ~= 3
    error('Input must have size 3 along the specified dimension.');
end

if (dim == 5)
    x0 = x;
    for c = 1:size(x, dim)
        for c_vol = 1:size(x, 4)
            x(:,:,:,c_vol,c) = medfilt3(x(:,:,:,c_vol,c), [5 5 5]);
        end
    end
end


% Sort the elements along the dimension
x_sorted = sort(x, dim);

% Get the three distances between pairs
x1 = slice_dim(x_sorted, dim, 1);
x2 = slice_dim(x_sorted, dim, 2);
x3 = slice_dim(x_sorted, dim, 3);

d12 = abs(x1 - x2);
d23 = abs(x2 - x3);
d13 = abs(x1 - x3);

% Compare distances
mask_12 = d12 <= d23 & d12 <= d13;
mask_23 = d23 < d12 & d23 <= d13;
mask_13 = d13 < d12 & d13 < d23;


% Sort the elements along the dimension
x_sorted = sort(x0, dim);

% Get the three distances between pairs
x1 = slice_dim(x_sorted, dim, 1);
x2 = slice_dim(x_sorted, dim, 2);
x3 = slice_dim(x_sorted, dim, 3);



% Compute averages where the mask applies
y = zeros_like(x1);
y(mask_12) = (x1(mask_12) + x2(mask_12)) / 2;
y(mask_23) = (x2(mask_23) + x3(mask_23)) / 2;
y(mask_13) = (x1(mask_13) + x3(mask_13)) / 2;
end

function s = slice_dim(x, dim, idx)
sz = size(x);
subs = repmat({':'}, 1, ndims(x));
subs{dim} = idx;
s = x(subs{:});
end

function z = zeros_like(x)
z = zeros(size(x));
end
