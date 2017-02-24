function LL = getLL_GRW_vectorized(params,sentData)

[ndtBounds,threshBounds,mu1IBounds,mu1SBounds,mu2IBounds,mu2SBounds,mu12IBounds,mu12SBounds] = getGRWParamBounds();

rtData = sentData(:,1);
choiceData = sentData(:,2);


% conditionInfo = sentData(:,3);
% isSpeed = conditionInfo > 6;
% isCued = mod(conditionInfo-1,6) >= 3;
% jitIndex = mod(conditionInfo,3) + 1;
% cueLevel = 

isSpeed = sentData(:,3);
isCued = sentData(:,4);
cueDev = sentData(:,5);
jitter = sentData(:,6);

mu1 = zeros(length(sentData),1);
mu2 = zeros(length(sentData),1);
ndt = zeros(length(sentData),1);
thresh = zeros(length(sentData),1);

mu1Intercept = infToDef(params(1),mu1IBounds(1),mu1IBounds(2));
mu1Slope = infToDef(params(2),mu1SBounds(1),mu1SBounds(2));
mu2Intercept = infToDef(params(3),mu2IBounds(1),mu2IBounds(2));
mu2Slope = infToDef(params(4),mu2SBounds(1),mu2SBounds(2));
ndt = infToDef(params(5),ndtBounds(1),ndtBounds(2));
threshSpeed = infToDef(params(6),threshBounds(1),threshBounds(2));
threshAcc = infToDef(params(7),threshBounds(1),threshBounds(2));
mu12Intercept = infToDef(params(8),mu12IBounds(1),mu12IBounds(2));
mu12Slope = infToDef(params(9),mu12SBounds(1),mu12SBounds(2));


muMagnitude = mu1Intercept + mu1Slope * jitter + mu12Intercept * isCued + mu12Slope * isCued .* abs(cueDev);
muAngle = isCued * mu2Intercept + mu2Slope * isCued .* cueDev;
mu1 = muMagnitude .* cos(muAngle);
mu2 = muMagnitude .* sin(muAngle);
thresh = isSpeed * threshSpeed + (1-isSpeed) * threshAcc;


sigmasq = 1;


k = 141; % sets precision of the infinite series used

j0k = besselzero(0,k,1);
J1 = besselj(1,j0k);

% the infinite series alternates around a number, cutting off at an even
% index will over estimate and cutting off at odd will under estimate.  So
% we can just divide the last term by two to get a rough approximation of
% the limit of the infinite series
J1(length(J1)) = J1(length(J1)) / 2;

%% Non-vectorized

rtL = 0;
for m = 1:length(j0k)
    rtL = rtL + (sigmasq) ./ (thresh.^2) .* ( (j0k(m) / J1(m)) .* exp(-1 * j0k(m)^2 * sigmasq * (rtData-ndt) ./ (2 * thresh.^2)));
end
% logrtLike = (rtData-ndt) * log(rtL);

rtL((rtData-ndt)<.01) = 10^(-100); % Get rid of weird model behavior that occurs when rt is too close to non-decision time


%logRespLike = (thresh./sigmasq) .* (mu1 .* cos(choiceData) + mu2 .* sin(choiceData)) - (1 / (2*sigmasq)) .* (mu1.^2 + mu2.^2) .* (rtData-ndt) ;


respLike = exp((thresh./sigmasq) .* (mu1 .* cos(choiceData) + mu2 .* sin(choiceData)) - (1 / (2*sigmasq)) .* (mu1.^2 + mu2.^2) .* (rtData-ndt)) ;

jointLikelihood = respLike .* rtL;

LL = -sum(logRespLike + log(rtL));

end
