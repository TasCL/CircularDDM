function [RT R] = rcircularddm2(n, threshold, angle, startpoint, t0, ...
                                p, tolerance)
%RCIRCULARDDM2  Generate random deviates for the circular
%drift-diffusion model
%
%   [RT R] = rcircularddm2(n, threshold, angle, startpoint, t0, p,
%   tolerance) calls rcddm2.mexa64 to generate random deviates for circular 
%   drift-diffusion distribution. Differing from rcircularddm1,
%   instead of taking a parameter vector with 5 numbers [a, vx, vy,
%   t0, s], rcircularddm2 takes a threshold vector, an angle
%   vector, a starting point matrix and a t0 scalar as model
%   parameters. Internally, in rcircularddm1, an angle is derived
%   from vx, vy and s. 
%
%   Input:
%     n          - number of observations.
%     threshold  - a user supplied threshold vector. It can be of
%                  length one or longer.
%     angle      - a user supplied angle vector. It can be of length one
%                  or longer
%     startpoint - a user suppled starting point matrix. First
%                  column is xPos and second column is yPos.
%     t0         - nondecision time. must be a scalar.
%     p          - a precision parameter for random walk step. Default
%                  is 0.15 seconds.
%     tolerance  - an upper bound for diffusion step. Default is
%                  1e3 steps.
%
%   Output:
%     RT   - random deviates for response times, 
%     R    - response angles.
%
%   Examples:
%     % Using an unrealistic example
%     t0 = 0.25;
%     threshold  = randn(1e3,1); 
%     angle      = rvonmises(1e3, 0, 10);
%     startpoint = rand(10,2)
%     tol        = 1e3;
%     [RT R] = rcircularddm2(20, threshold, angle, startpoint, t0, 0.15, tol);
%     [RT(1:10,:) R(1:10,:)]  % Show the first 10 rows
% 
%     figure(1)
%     histogram(RT)
%     xlabel('Response time')
%     
%     figure(2)
%     histogram(R)
%     xlabel('Responses') 
%     
%     help('rcircularddm2') % Show this help page
%
%   References:
%     Smith, P. L. (2016). Diffusion Theory of Decision Making in
%        Continuous Report, Psychological Review, 123 (4),
%        425--451.
%
% Circular Drift-diffusion Model for Unix-like MATLAB
% (c) Yi-Shin Lin, 2017, yishin.lin@utas.edu.au
if nargin < 6
     p = 0.15;
     tolerance = 1e3;
end

[RT R] = rcddm2(n, threshold, angle, startpoint, t0, p, tolerance);
