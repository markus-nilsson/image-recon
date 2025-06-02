function O = imop_harmonize(P, T, do_use_background_from_target)

if (nargin < 3), do_use_background_from_target = 1; end

% Address bad predictions
P(P < 0) = 0;

M = median(P,5);

% Mask it
for c_vol = 1:size(M,4)
    tmp1 = mio_mask_threshold(M(:,:,:,c_vol));
    tmp2 = mio_mask_threshold(T(:,:,:,c_vol));

    tmp = tmp1 & tmp2;

    tmp = mio_mask_keep_largest(tmp);
    
    mask(:,:,:,c_vol) = tmp;
end

mask = mean(double(mask),4) > 0.8;
mask = mio_mask_erode(mask);
mask = mio_mask_keep_largest(mask);
mask = mio_mask_fill(mask);


% rescale
for c_vol = 1:size(M,4)

    A = M(:,:,:,c_vol);
    B = T(:,:,:,c_vol);

    
    % rough alignment first
    A = A / (B(mask(:)) \ A(mask(:)));

    % refine alignment
    R = A(mask(:)) ./ B(mask(:));
    s = median(R( (R > 0.3) & (R < 2)));
    
    A = A / s;

    % fit a polynomial ratio
    R = fitRatioVolume(B .* mask, A .* mask, 'Polyorder', 3);

    R = mio_min_max_cut(R, [0.1 10]);

    % Multipy 
    A = A .* R;

    % Use background from target
    if  (do_use_background_from_target)
        A(mask(:) == 0) = B(mask(:) == 0);
    end

    M(:,:,:,c_vol) = A;

end

O = M;
