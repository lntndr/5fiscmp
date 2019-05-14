function [s,sstr]=lorenzsign(rho,N)
%LORENZSIGN Computes an N-elements signature of the trajectory in a
% Lorenz attractor starting from yc0 + [0; 0; 3], where
% yc0 = [rho-1; eta; eta], eta = sqrt(beta*(rho-1)), sigma=10, beta=8/3 and
% the given rho. In this set of hypothesis every trajectory reaches a
% series of local maxima in y(1). The signature is created by associating
% to all maxima '0' if y(2:3)>0, '1' otherwise.
%
%          S = LORENZSIGN(RHO,N) returns in s a N-row long
%   array that represents the signature.
%
%   [S,SSTR] = LORENZSIGN(RHO,N) returns in s a logical N-row long
%   array and in sstr the same array split in substrings of fixed lenght to
%   improve readibility truncating the last k element, k=mod(N,sstrl).
%
%   Expected runtime r with console output suppressed by ;
%   if    rho<470/19      r < 1e-3 s
%   else                  r < 4e-3*N s
%   (Extimated with MATLAB 2019a, sum(bench)=2.395+-0.008)

narginchk(2,2);

sstrl = 60;

if rho < 1
    error('ErrorTAG:TagName', strcat ( 'Signature is not defined if', ...
        ' rho < 1 as [0 0 0] is the only critical point.') );
elseif rho == 1
    error('ErrorTAG:TagName', strcat ( 'The function is not designed', ...
        ' to deal with the simple bifurcation for rho = 1.') );
elseif rho < 470/19
    % If 1 < rho <470/19 yc0 and yc1 are globally stable so for the given
    % y0 the orbit never stops to wind around yc0.
    s=zeros(N,1);
else
    % Display warning for huge rho values
    if rho > 1e7
        % For rho > 1e7 the part of the function that enumerates events
        % unexpected behaviour, namely it stops the ode solver too late.
        % The function still shows the pseudo-periodic behaviour expected
        % for large values of rho, but it can not guarantee the
        % correctness of results.
        warning('The function is not designed to deal with rho > 1e7.');
    end
    
    % Initialize system parameters
    sgm = 10;
    bta = 8/3;
    eta = sqrt(bta*(rho-1));
    A = [-bta, 0, eta; 0, -sgm, sgm; -eta, rho, -1];
    yc = [rho-1; eta; eta];
    y0 = yc + [0; 0; 3];
    tme = [0 Inf]; % Check localy1max to see when ode113 stops

    % Initialize parameters for the nested event function localy1max
    k=(N-1)*2; 
        %The method localy1max uses is unable to distinguish local minima
        %from local maxima, so it looks for the double of the events.
    pd=0; %Previous difference
    nd=0; %New difference
    py1=0; %Previous y(1)
    j=0;

    opts = odeset('Events',@localy1max,'RelTol',1.e-10,'AbsTol',1.e-15);
    [~,~,~,ye]=ode113(@lorenzeqn, tme, y0, opts, A);
    
    s=(ye(:,2)<0);
    
end

if nargout>1 %sstr is parsed only if needed
    sstr=regexp(sprintf('%d', s), sprintf('\\w{1,%d}', sstrl), 'match')';
    if size(sstr{end},2)<sstrl
        sstr(end)=[];
    end
end

    function [value,isterminal,direction] = localy1max(~,y,A)
    % LOCALY1MAX stops the ODE solver after the solution has reached N 
    % local maxima. As MATLAB does not allow to check out the size of ye 
    % array, this function works around defining two events: the first one 
    % is a zero cross of the derivative, calculated internally by odeevents
    % , the second one is a counter reduced by one every time the function
    % reaches a local critical point. 
    % This function is intentend to be nested into lorenzsign.
    
    ydot = lorenzeqn(y,y,A); %The first entry is only a placeholder
    
    if mod(j,5)==0
        % 5 is a sampling value obtained experimentally, supposed to be
        % good in the range 470/19 < rho < 1e7. Sampling is aimed to avoid
        % false positives caused by numerical noise.
        
        nd=y(1)-py1; % Evaluates the difference between the the current 
                     % y(1) and the previous one
        if pd*nd<0
            % If the sign changes beetween the current and the previous
            % difference, then a local critical point is reached
            k=k-1;
        end
        % Saves the current values for the next cycle
        pd=nd;
        py1=y(1);
    end
    
    j=j+1;
    
    value = [ydot(1) k];
    isterminal = [0 1];
    direction = [-1 0];
    
    end

end

function ydot = lorenzeqn(~,y,A)
%LORENZEQN  Equation of the Lorenz chaotic attractor.
%   ydot = lorenzeqn(t,y,A).
%   The differential equation is written in almost linear form.
%      ydot = A*y
%   where
%      A = [ -beta    0     y(2)
%               0  -sigma   sigma 
%            -y(2)   rho    -1  ];

A(1,3) = y(2);
A(3,1) = -y(2);
ydot = A*y;
end