
function E_B = NAA_optimizer(Nrisk, mpc, N_Bess, P_profile, P_renewable)
    rng(123);
    % Algorithm parameters
    pop_size = 20;      % number of individuals in the population
    max_iter = 100;     % maximum number of iterations (generations)
    E_min = 0.5;        % lower bound for capacity (MWh)
    E_max = 5;          % upper bound for capacity (MWh)
    %penalty_weight = 1e4; % weight for voltage violation penalty in the cost


    % 1. Initialize Population
    % Each row in 'pop' corresponds to a candidate vector of BESS capacities.
    pop = E_min + (E_max - E_min) .* rand(pop_size, N_Bess);
    fitness_array = zeros(pop_size, 1);
    
    for i = 1:pop_size

        E_B = pop(i,:);
        cost = TotalCost(E_B, mpc, Nrisk, P_profile, P_renewable);
        fitness_array(i) = cost;
    end
    
    % 2. Shelter Initialization
    numShelters = 5;
    [shelters, shelterLeaders] = shelterInitialization(fitness_array, numShelters);

    E_B = pop(1, :);

    assignin("base","pop", pop);
    assignin("base","shelters",shelters);
    assignin("base","shelterLeaders",shelterLeaders);
    

end

function [shelters] = shelterInitialization(fitness, numShelters)
    
    [~, sortedIndices] = sort(fitness, 'ascend');  % Sorted indices so that lower fitness is better
    
    % Initialize a cell array to hold candidate indices for each shelter
    shelters = cell(numShelters, 1);
    for i = 1:numShelters
        shelters{i} = [];  % Initialize each shelter as an empty array [2]
    end
    
    % Assign candidates to shelters in round-robin fashion using the sorted indices
    for idx = 1:length(sortedIndices)
        shelterIdx = mod(idx-1, numShelters) + 1;  % Determine the shelter index (round-robin assignment)
        shelters{shelterIdx} = [shelters{shelterIdx}, sortedIndices(idx)];
    end
    
    % In each shelter, find the candidate with the lowest fitness value to be the leader
    shelterLeaders = zeros(numShelters, 1);
    for i = 1:numShelters
        currentShelter = shelters{i};
        [~, bestLocalIdx] = min(fitness(currentShelter));  % Identify the best candidate in the shelter [2]
        shelterLeaders(i) = currentShelter(bestLocalIdx);
    end
end
