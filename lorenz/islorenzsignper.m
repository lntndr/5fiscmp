function per=islorenzsignper(rho,N,lastn)

narginchk(1,3);

if nargin < 3
    lastn=50;
end

if nargin < 2
    N=200;
end

[a,~] = lorenzsign(rho,N);

if isempty(a(find(a,1):end))
    disp("With the given parameters the trajectory revolves around y0");
    per=1;
else
    a = a(end-lastn:end); %remove transient
    xa = xcorr(a);
    xa = xa(ceil(length(a)/2):length(a));
    [~,i] = findpeaks(xa);
    pv = diff(i);
    if ~range(pv)
        per = pv(1);
    else
        per = inf;
    end
end