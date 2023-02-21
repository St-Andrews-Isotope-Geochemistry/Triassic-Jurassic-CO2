% Reads in evolutions and shapes them into a structure of arrays
function output = getShapedEvolutions(filename)
    raw_evolutions = readmatrix(filename);
    reshaped_evolutions = reshape(raw_evolutions,[100,13,numel(raw_evolutions)/(13*100)]);

    output.pH = squeeze(reshaped_evolutions(:,1,:));
    output.co2 = squeeze(reshaped_evolutions(:,2,:))*1e6;
    output.saturation_state = squeeze(reshaped_evolutions(:,3,:));
    output.dic = squeeze(reshaped_evolutions(:,4,:));
    output.alkalinity = squeeze(reshaped_evolutions(:,5,:));
    output.temperature = squeeze(reshaped_evolutions(:,6,:));
    output.d11B = squeeze(reshaped_evolutions(:,7,:));
    output.calcium = squeeze(reshaped_evolutions(:,8,:));
    output.magnesium = squeeze(reshaped_evolutions(:,9,:));
    output.epsilon = squeeze(reshaped_evolutions(:,10,:));
    output.d11B_sw = squeeze(reshaped_evolutions(:,1,:));
    
end