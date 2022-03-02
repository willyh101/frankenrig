function out = idx2dirs(visIDs)

% assumes that 1 is a blank trial
ndirs = max(visIDs);
dirs = linspace(0,360,ndirs);
dirs = [nan dirs(1:end-1)];
out = arrayfun(@(x) dirs(x), visIDs);