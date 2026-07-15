function spikedata = fixSpikes(spikedata)
spikedata = spikedata{1};
for j=1:length(spikedata)
    if isequal(spikedata{j},-1)
        spikedata{j} = [];
    else
        sd = struct(spikedata{j});
        spikedata{j} = sd.t/10000;
    end
end