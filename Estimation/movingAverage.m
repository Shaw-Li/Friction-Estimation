function mvAvg = movingAverage(input, window)
%movingAverage: This function calculates the moving average for the given
%data

mvAvg = zeros(length(input) - window, 1);

for i=1:length(mvAvg)
    mvAvg(i) = mean(input(i:i+window));
end

end

