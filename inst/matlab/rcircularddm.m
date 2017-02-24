function [RT R A] = rcircularddm(n, pVec, p)
%RCIRCULARDDM  Generate random deviates for the circular
%drift-diffusion model
%
%   [RT R A] = rcircularddm(n, pVec, stepTime) calls rddm.mexa64 to
%   generate random deviates for circular drift-diffusion
%   distribution.
%
%   Input:
%     n    - number of observations.
%     pVec - a circulr DDM parameter vector with the order, a, vx,
%            vy, t0, and s. Each represents decision threshold,
%            drift rate for x axis, drift rate for y axis,
%            nondecision time and scaling parameter.
%     p    - a precision parameter for random walk step. Default
%             is 0.15 seconds.
%
%   Output:
%     RT   - random deviates for response times, 
%     R    - continuous reports
%     A    - response angles.
%
%   Examples:
%     % threshold=2; vx=1.5; vy=1.25; t0=0.25; sigma_square = 1;
%     pVec     = [2, 1.5, 1.25, .25, 1]; 
%     stepTime = .001;  % use 1 ms step time, instead of 0.15 s
%     [RT R A] = rcircularddm(1e3, pVec, stepTime);
%     
%     [RT(1:10,:) R(1:10,:) A(1:10,:) ]  % Show the first 10 rows
% 
%     figure(3)
%     histogram(RT)
%     xlabel('Response time')
%     
%     figure(4)
%     histogram(A)
%     xlabel('Response angle')
%     
%     help('rcircularddm') % Show this help page
%
%   References:
%     Smith, P. L. (2016). Diffusion Theory of Decision Making in
%        Continuous Report, Psychological Review, 123 (4),
%        425--451.
%
% Circular Drift-diffusion Model for Unix-like MATLAB
% (c) Yi-Shin Lin, 2017, yishin.lin@utas.edu.au
if nargin < 3
     p = 0.15;
end

[RT R A] = rcddm(n, pVec, p);
