% encode string information
data.ori = 300;
data.contrast = 100.3;
data.position = [0 100];
data.size = 50;
%%
s = jsonencode(data);
temp = reshape(dec2hex(s),1,[]);
out = hexToBinaryVector(temp);

%% check it
invar = binaryVectorToHex(out);
recvd = char(hex2dec(reshape(invar, [], 2)))';
the_data = jsondecode(recvd);