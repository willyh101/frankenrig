function result = decodeJSONbin(signal)

% get hex
hexStr = binaryVectorToHex(signal);

% hex to json string
json_str = char(hex2dec(reshape(hexStr, [], 2)))';
result = jsondecode(json_str);