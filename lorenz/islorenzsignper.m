function per=islorenzsignper(rho,N,lastn)
% ISLORENZSIGNPER determines if the signature of the trajectory of the
% Lorenz attactor associated to rho is periodic. Depending on the input, it
% can compute the signature or merely look for the period.

narginchk(1,3);

% Default settings
if nargin < 3
    lastn=50;
end

if nargin < 2
    N=200;
end

if length(rho) == 1 % Density associated to a signature
    [a,~] = lorenzsign(rho,N);
else % Already a signature
    a = rho;
end

if isempty(a(find(a,1):end))
    disp("With the given parameters the trajectory revolves around y0");
    per=1;
else
    a = a(end-lastn:end); %remove transient
    xa = xcorr(a); %Autocorrelation of the signature
    xa = xa(ceil(length(a)/2):length(a));
    [~,i] = findpeaks(xa); % Maxima of correlation are hints of periodicity
    pv = diff(i); % Find the distance between maxima
    if ~range(pv) % If all the maxima are equidistant, there's a period
        per = pv(1);
    else
        per = inf;
    end
end