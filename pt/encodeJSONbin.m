function out = encodeJSONbin(data)

% encode into JSON string
s = jsonencode(data);

% make into hex and then move to binary
temp = reshape(dec2hex(s),1,[]);
out = hexToBinaryVector(temp);