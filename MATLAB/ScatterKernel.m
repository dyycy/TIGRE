function [outputArg1,outputArg2] = ScatterKernel(inputArg1,inputArg2)
%SCATTERKERNEL Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end



function ce = kernel_ce(Prm, Blk, du, dv, A, alpha, beta)

edgewt = imbinarize(Prm,T);
% ce = A .* edgewt .* (Prm./Blk).^(alpha) .* (log(Blk./(Prm+eps))).^(beta);

