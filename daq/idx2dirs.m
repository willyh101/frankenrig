function out = idx2dirs(visIDs, has_blank)

if nargin < 2
    has_blank = 1;
end

% assumes that 1 is a blank trial
ndirs = max(visIDs);
dirs = linspace(0,360,ndirs);
if has_blank
    dirs = [nan dirs(1:end-1)];
else
    dirs = dirs(1:end-1);
end
out = arrayfun(@(x) dirs(x), visIDs);