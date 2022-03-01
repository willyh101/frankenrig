repNum = 50;
repFields = [7]; [10 11 12];
nonRepNumbers =repFields(end)+1:repFields(end)+repNum;

list = [repFields nonRepNumbers repmat(repFields,[1 repNum-1])];

%%
