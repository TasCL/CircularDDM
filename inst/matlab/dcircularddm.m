function out = dcircularddm(x, pVec, k)
%DCIRCULARDDM  Log-likelihood function for the circular
%drift-diffusion model
%   out = dcircularddm(x, pVec, k) calls dcddm.mexa64 to calculate
%   log-likelihood for circular drift-diffusion distribution.
%
%   Input:
%     x     - a matrix with a response-time column (1st) and a
%     continuous outcome column. Each row represents a trial.
%     pVec  - a circulr DDM parameter vector with the order, a, vx,
%     vy, t0, and s. Each represents decision threshold, drift rate
%     for x axis, drift rate for y axis, nondecision time and
%     scaling parameter.
%     k     - a precision parameter for calculating the infinite
%     series. The larger the k, the more memory space is required
%     to conduct the computation. Default is 141.
%
%   Output:
%     out   - random deviates from von Mises distribution
%
%   Examples:
%     X = [1.2595272 1.9217430; % A 6-by-2 matrix
%          0.8693937 1.7844653; 
%          0.8009044 0.2662521; 
%          1.0018933 2.1569724;
%          2.3640007 1.7277440;
%          1.0521304 0.8607271];
%     pVec = [2.45, 1.5, 1.25, .1, 1]; 
%     dcircularddm(X, pVec) % LL fucntion for circular DDM
%     help('dcircularddm')
%
%   References:
%     Smith, P. L. (2016). Diffusion Theory of Decision Making in
%        Continuous Report, Psychological Review, 123 (4),
%        425--451.
%
% Circular Drift-diffusion Model for Unix-like MATLAB
% (c) Yi-Shin Lin, 2017, yishin.lin@utas.edu.au
if nargin < 3
     k = 141;
end

out = dcddm(x, pVec, k);
