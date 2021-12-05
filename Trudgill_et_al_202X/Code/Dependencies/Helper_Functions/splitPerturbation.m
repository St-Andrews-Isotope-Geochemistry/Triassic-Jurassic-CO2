% Function to split preperturbation from during perturbation
function [preperturbation,perturbation] = splitPerturbation(input,boundary)
    preperturbation_boolean = input.age>=boundary;
    perturbation_boolean = input.age<boundary;
    
    preperturbation.pH = reshape(input.pH(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.co2 = reshape(input.co2(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.saturation_state = reshape(input.saturation_state(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.dic = reshape(input.dic(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.alkalinity = reshape(input.alkalinity(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.temperature = reshape(input.temperature(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.d11B = reshape(input.d11B(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.calcium = reshape(input.calcium(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.magnesium = reshape(input.magnesium(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.epsilon = reshape(input.epsilon(preperturbation_boolean),[],size(input.pH,2));
    preperturbation.d11B_sw = reshape(input.d11B_sw(preperturbation_boolean),[],size(input.pH,2));

    perturbation.pH = reshape(input.pH(perturbation_boolean),[],size(input.pH,2));
    perturbation.co2 = reshape(input.co2(perturbation_boolean),[],size(input.pH,2));
    perturbation.saturation_state = reshape(input.saturation_state(perturbation_boolean),[],size(input.pH,2));
    perturbation.dic = reshape(input.dic(perturbation_boolean),[],size(input.pH,2));
    perturbation.alkalinity = reshape(input.alkalinity(perturbation_boolean),[],size(input.pH,2));
    perturbation.temperature = reshape(input.temperature(perturbation_boolean),[],size(input.pH,2));
    perturbation.d11B = reshape(input.d11B(perturbation_boolean),[],size(input.pH,2));
    perturbation.calcium = reshape(input.calcium(perturbation_boolean),[],size(input.pH,2));
    perturbation.magnesium = reshape(input.magnesium(perturbation_boolean),[],size(input.pH,2));
    perturbation.epsilon = reshape(input.epsilon(perturbation_boolean),[],size(input.pH,2));
    perturbation.d11B_sw = reshape(input.d11B_sw(perturbation_boolean),[],size(input.pH,2));
end