function val = decodeBinaryDaqSignal(clockSweep, pulseSweep)

clk_times = find(diff(clockSweep) == 1);
stm_times = find(diff(pulseSweep) == 1);
binvec = ismember(clk_times, stm_times)';
val = binaryVectorToDecimal(binvec);