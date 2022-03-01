function [w h] = splitPretty(N,L,H,verbose);
%split N into integers that fit the ratio of L to H
% if L<H %ensure long is always longer
%     temp =H;
%     H=L;
%     L=temp;
% end

if nargin<4
    verbose =1;
end

b1 = 0.5;
b2 = 2;
b3 = 1;

fun = @(x) b1*sum(abs(x)) + b2*abs(x(1)/x(2) - L/H) + b3*prod(x);
% fun = @(x) b2*abs(x(1)/x(2) - L/H) + b3*prod(x);

% A = [1 1];
% b = N;

lb = [0 0];
ub = [N N];

% c = @(x) deal(1000*(N-x(1)*x(2)),[]); %don't allow the product of the two to be less than the desired N
% c = @(x) deal([],floor((N-0.01)/x(1)*x(2))); %don't allow the product of the two to be less than the desired N
c = @(x) deal(floor((N-0.01)/(x(1)*x(2)))+(N-x(1)*x(2)),[]); %don't allow the product of the two to be less than the desired N

options = optimoptions('ga');
options.FunctionTolerance = 1e-6; %default 1e-6
options.MaxStallGenerations = 50;
if verbose==0
    options.Display='off';
end
% [x fVal] = ga(fun,2,A,b,[],[],lb,ub,c,[1 2])
[x fVal] = ga(fun,2,[],[],[],[],lb,ub,c,[1 2],options);

x = max(x,[1 1]);


w = x(2);
h = x(1);

if verbose
disp(['WxH: ' num2str(w) 'x' num2str(h) '. Ratio: ' num2str(w/h,3) '. Ideal: ' num2str(L/H,3) '. Product: ' num2str(w*h)]);
end
