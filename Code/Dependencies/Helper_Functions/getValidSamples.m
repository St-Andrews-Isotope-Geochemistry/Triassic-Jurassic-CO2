%% Function that returns only valid samples
function output = getValidSamples(input)
    time_series_length = size(input.co2,1);
    boolean = ~(repmat(any(isnan(input.pH)),time_series_length,1) | repmat(input.co2(end,:)<500,time_series_length,1) | repmat(input.co2(end,:)>5000,time_series_length,1));
    output.pH = reshape(input.pH(boolean),time_series_length,[]);
    output.co2 = reshape(input.co2(boolean),time_series_length,[]);
    output.saturation_state = reshape(input.saturation_state(boolean),time_series_length,[]);
    output.dic = reshape(input.dic(boolean),time_series_length,[]);
    output.alkalinity = reshape(input.alkalinity(boolean),time_series_length,[]);
    output.temperature = reshape(input.temperature(boolean),time_series_length,[]);
    output.d11B = reshape(input.d11B(boolean),time_series_length,[]);
    output.calcium = reshape(input.calcium(boolean),time_series_length,[]);
    output.magnesium = reshape(input.magnesium(boolean),time_series_length,[]);
    output.epsilon = reshape(input.epsilon(boolean),time_series_length,[]);
    output.d11B_sw = reshape(input.d11B_sw(boolean),time_series_length,[]);
end