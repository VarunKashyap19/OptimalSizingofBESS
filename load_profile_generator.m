function [p_profile,q_profile] = load_profile_generator(mpc_bus)
    rng(123);
    %baseMultipliers = [0.40, 0.38, 0.35, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.90, 0.85, 0.80, 0.85, 0.90, 0.96, 0.96, 0.95, 0.85, 0.75, 0.65, 0.50];
    baseMultipliers = [0.70, 0.68, 0.65, 0.65, 0.70, 0.90, 1.20, 1.40, 1.60, 1.70, 1.80, 1.90, 1.80, 1.70, 1.60, 1.70, 1.80, 1.92, 1.92, 1.90, 1.70, 1.50, 1.30, 1.00];
    variation_range = 0.04;

    num_buses = size(mpc_bus, 1);
    base_P = mpc_bus(:, 3);
    base_Q = mpc_bus(:, 4); %shayad
    num_hours = 24;

    p_profile = zeros((num_buses), 24);
    q_profile = zeros((num_buses), 24);

    randomizedMultipliers = zeros(num_buses, num_hours);
    q_randomizedMultipliers = zeros(num_buses, num_hours);

    % Randomized multipliers generation for P
    for bus = 1:num_buses
        variations = (rand(1, num_hours) * 2 - 1) * variation_range;
        randomizedMultipliers(bus, :) = baseMultipliers + variations;
    end

    % Generate P profile for each bus
    for bus = 1:num_buses
        p_profile(bus, :) = base_P(bus) * randomizedMultipliers(bus, :);
    end
    

    q_baseMultipliers = [0.39, 0.36, 0.34, 0.34, 0.39, 0.49, 0.59, 0.72, 0.82, 0.87, 0.92, 0.97, 0.92, 0.87, 0.82, 0.87, 0.92, 1.00, 1.00, 0.97, 0.87, 0.74, 0.64, 0.49];
    
    % Randomized multipliers generation for Q
    for bus = 1:num_buses
        variations = (rand(1, num_hours) * 2 - 1) * variation_range;
        q_randomizedMultipliers(bus, :) = q_baseMultipliers + variations;
    end


    % Generate Q profile for each bus
    for bus = 1:num_buses
        q_profile(bus, :) = base_Q(bus) * q_randomizedMultipliers(bus, :);
    end
    
end