function [fittedVolume, maskFinal] = imop_fit_ratio(vol1, vol2, varargin)
%FITRATIOVOLUME Fits a smooth 3D polynomial to the ratio of vol1 ./ vol2
%   [fittedVolume, maskFinal] = fitRatioVolume(vol1, vol2, ...)
%
%   Inputs:
%       vol1, vol2       - 3D arrays of same size
%
%   Optional parameters:
%       'PolyOrder'      - Order of polynomial fit (default: 2)
%       'MaxIter'        - Max number of iterations (default: 5)
%       'SigmaThresh'    - Outlier threshold in standard deviations (default: 2)
%
%   Outputs:
%       fittedVolume     - Fitted polynomial values
%       maskFinal        - Final inlier mask

% Parse input
p = inputParser;
addParameter(p, 'PolyOrder', 2);
addParameter(p, 'MaxIter', 5);
addParameter(p, 'SigmaThresh', 2);
parse(p, varargin{:});

order = p.Results.PolyOrder;
maxIter = p.Results.MaxIter;
sigmaThresh = p.Results.SigmaThresh;

% Create ratio

vol1(vol1 < 0) = 0;
vol2(vol2 < 0) = 0;

ind_bg = (~mio_mask_threshold(vol1)) | (~mio_mask_threshold(vol2));

ind_bg = ind_bg | (vol1 < eps) | (vol2 < eps);

tol1 = vol1; tol1(ind_bg) = quantile(tol1(~ind_bg), 0.01);
tol2 = vol2; tol2(ind_bg) = quantile(tol2(~ind_bg), 0.01);


ratio = tol1 ./ tol2;
ratio = medfilt3(ratio, [3 3 3]);

mask = (~ind_bg) & isfinite(ratio) & (vol2 ~= 0) & (vol1 ~= 0);

% Eliminate really off ratios
r_upper = median(ratio(mask(:))) + 5 * mad(ratio(mask(:)));
r_lower = median(ratio(mask(:))) - 5 * mad(ratio(mask(:)));
mask = mask & (ratio < r_upper) & (ratio > r_lower);


sz = size(vol1);
[x, y, z] = ndgrid(1:size(vol1,1), 1:size(vol1,2), 1:size(vol1,3));
x = x(:) / sz(1); y = y(:) / sz(2); z = z(:) / sz(3);  % normalize to [0,1]
x = x(:); y = y(:); z = z(:); ratioVec = ratio(:); maskVec = mask(:);
X0 = [x(maskVec), y(maskVec), z(maskVec)];
Y0 = ratioVec(maskVec);


% Construct polynomial terms
polyTerms = @(X) buildPolyTerms(X, order);

% Iterative fitting
X = X0;
Y = Y0; 
for iter = 1:maxIter
    A = polyTerms(X);
    coeffs = A \ Y;
    Yfit = A * coeffs;
    residuals = Y - Yfit;
    stdRes = mad(residuals) / 0.8;
    keep = abs(residuals) < sigmaThresh * stdRes;
    X = X0(keep, :);
    Y = Y0(keep);
end

% Evaluate on full grid
Xfull = [x(:), y(:), z(:)];
Afull = polyTerms(Xfull);
Yfull = reshape(Afull * coeffs, size(vol1));
fittedVolume = Yfull;
maskFinal = false(size(vol1));
maskFinal(mask) = true;
% maskFinal(mask) = keep;

end

function A = buildPolyTerms(X, order)
% Build polynomial design matrix up to specified order
% E.g. for order = 2: 1, x, y, z, x^2, xy, xz, y^2, yz, z^2
x = X(:,1); y = X(:,2); z = X(:,3);
A = ones(size(x));
for i = 0:order
    for j = 0:(order - i)
        for k = 0:(order - i - j)
            if i + j + k > 0
                A = [A, (x.^i).*(y.^j).*(z.^k)];
            end
        end
    end
end
end
