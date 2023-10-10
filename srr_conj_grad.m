function [x,cgr] = srr_conj_grad(E,b,tol,maxit,x)
% A conjugate gradient solver that allows to plot the result in each
% iteration
% Solves: E*x=b
%
% x = conjugate_gradient(E,b)
% x = conjugate_gradient(E,b,tol)
% x = conjugate_gradient(E,b,tol,maxit)
% x = conjugate_gradient(E,b,tol,maxit,z)
% x = conjugate_gradient(E,b,tol,maxit,z,verbose)
% x = conjugate_gradient(E,b,tol,maxit,z,verbose,newfigure)
% [x,x_iter] = conjugate_gradient(_____)
%
%
% Input:
%   E         =  Hermitian and positive definite matrix or operator that
%                immitates such a matrix and implements a mtimes funciton.
%                Alternatively, A can function handle that imitates such a
%                matrix.
%   b         =  Right hand side of the equation A*x=b
%   tol       =  Tolerance of the method (default: 1e-6)
%   maxit     =  maximum number of iterations (default: 20)
%   x         =  initial guess (default: 0)
%   verbose   =  0 for no output, 1 for plotting the images and in each
%                iteration and 2 for printing the residal of the for each
%                iteration.
%   newfigure =  0 for using the same figure for each call of the funciton
%                or 1 for opening a new figure.
%
%
% Output:
%   x      = Result
%   x_iter = Result after each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col = @(x) x(:);
% norm = @(x) sum(col(x).^2);

if (1)
    Ed = E(x);
    d = b - Ed;
    r = d;
    
    normrr0 = norm(b);%.norm();
    normrr  = norm(r);%r.norm();
    
    for n = 1:maxit
        
        Ed = E(d);
        tmp = d .* Ed;


        alpha = normrr/(sum(col(tmp))+eps);
        x = x + alpha*d;

        r = r - alpha*Ed;
        normrr2 = norm(r) + eps;%norm();
        beta = normrr2/normrr;
        normrr = normrr2;
        d = r + beta*d;

        if (1)
            tmp = x.imreshape();
            subplot(1,3,1);  msf_imagesc(tmp(:,:,:,1));
            subplot(1,3,2);  msf_imagesc(tmp(:,:,:,2));
            subplot(1,3,3);  msf_imagesc(tmp(:,:,:,3));
            pause(0.05);

        end

        cgr(n) = sqrt(normrr/normrr0);

        if (cgr(n) < tol), break; end
        
    end

end
