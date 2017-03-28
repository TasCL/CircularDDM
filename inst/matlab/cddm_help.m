help('rcircularddm')
help('dcircularddm')
help('rvonmises')

clear all; cd('/home/yslin/Documents/MATLAB/');
mex rcddm.cpp -I/usr/include/gsl/ -larmadillo -lgsl

   % threshold=2; vx=1.5; vy=1.25; t0=0.25; sigma_square = 1;
      pVec     = [2, 1.5, 1.25, .25, 1]; 
      stepTime = .001;  % use 1 ms step time, instead of 0.15 s
      [RT R A] = rcircularddm(1e3, pVec, stepTime);
      
      [RT(1:10,:) R(1:10,:) A(1:10,:) ]  % Show the first 10 rows
  
      figure(3)
      histogram(RT)
      xlabel('Response time')
      
      figure(4)
      histogram(A)
      xlabel('Response angle')