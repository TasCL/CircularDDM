function out = rvonmises(n, mu, kappa)
%RVONMISES  Random number generator for the von Mises circular
%distribution.
%   out = rvonmises(n, mu, kappa) calls rvm.mexa64 to generates
%   random deviates for a von Mises distribution.
%
%   Input:
%     n     - number of observations 
%     mu    - mean direction of the distribution
%     kappa - width
%
%   Output:
%     out   - random deviates from a von Mises distribution
%
%   Examples: 
%     n       = 20; 
%     mu      = 0; 
%     kappa   = 10;
%     obs_deg = rvonmises(n, mu, kappa);
%     obs_rad = mod(obs_deg, 2*pi);
%
%   References:
%     Fisher, N. I., (1993). Statistical analysis of circular
%        data, Section 3.3.6, p. 49.
%     Jammalamadaka, S. R. & SenGupta, A. (2001). Topics in
%        Circular Statistics, Section 2.2.4, World Scientific
%        Press, Singapore.
%
% Circular Drift-diffusion Model for Unix-like MATLAB
% (c) Yi-Shin Lin, 2017, yishin.lin@utas.edu.au
out = rvm(n, mu, kappa);
