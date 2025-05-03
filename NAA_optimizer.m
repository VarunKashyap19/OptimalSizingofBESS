
function E_B = NAA_optimizer(Nrisk, mpc, N_Bess, P_profile, P_renewable)
    rng(123);
    % Algorithm parameters
    pop_size = 20;      % number of individuals in the population
    max_iter = 20;     % maximum number of iterations (generations)
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

    % 3. Stochastic Migration Model
    migrationRate = 0.8;

    fprintf('      Initial population Results.\n');
    fprintf('      Fitness:\n');
    disp(fitness_array);
    fprintf('      Leaders:\n');
    disp(shelterLeaders);

    for i = 1:max_iter
        [pop, fitness_array_temp, shelters, shelterLeaders_temp] = stochasticMigration(pop, fitness_array, shelters, shelterLeaders, mpc, Nrisk, P_profile, P_renewable, E_min, E_max, migrationRate);
        fprintf('\n\n      ******************************************\n');
        fprintf('      ******************************************\n');
        fprintf('      Changes after %d generation.\n', i);
        fprintf('      Cost Changes:\n');
        change = fitness_array - fitness_array_temp;
        disp(change);
        fprintf('      New Leaders:\n');
        disp(shelterLeaders_temp);
        fprintf('      ******************************************\n');
        fprintf('      ******************************************\n\n\n');
        fitness_array = fitness_array_temp;
        shelterLeaders = shelterLeaders_temp;
    end





    E_B = pop(1, :);

    assignin("base","pop", pop);
    assignin("base","fitness_array",fitness_array);
    assignin("base","shelters",shelters);
    assignin("base","shelterLeaders",shelterLeaders);
    

end

function [shelters, shelterLeaders] = shelterInitialization(fitness, numShelters)
    
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

    shelterLeaders = sortedIndices(1:numShelters);
end


function [pop, fitness_array, shelters, shelterLeaders] = stochasticMigration(pop, fitness_array, shelters, shelterLeaders, mpc, Nrisk, P_profile, P_renewable, E_min, E_max, migrationRate)
    pop_size = size(pop, 1);
    numVars = size(pop, 2);
    numShelters = length(shelters);

    max_deviation = fitness_array(shelters{end}(end)) - fitness_array(shelters{1}(1));

    for candidate = 1:pop_size
        currentShelterIdx = find(cellfun(@(s) any(s == candidate), shelters));

        % should not occur if every candidate is assigned.

        currentShelter = currentShelterIdx(1);
        currLeaderIdx = shelterLeaders(currentShelter);
        leaderFitness = fitness_array(currLeaderIdx);

        if fitness_array(candidate) > leaderFitness
            L_p = migrationRate * ((fitness_array(candidate) - leaderFitness) / max_deviation);
        else
            L_p = 0;
        end

%        fprintf('\nCandidate %d in shelter %d: Fitness = %.4f, Leader Fitness = %.4f, L_p = %.4f\n',candidate, currentShelter, fitness_array(candidate), leaderFitness, L_p);

        if rand < L_p
            % Explorer Branch
%            fprintf('  -> Candidate %d selected for EXPLORER branch.\n', candidate);
            pos = find(shelters{currentShelter} == candidate);
            shelters{currentShelter}(pos) = [];
        
            Fitness = fitness_array(candidate);
        
            while true
                % Randomly select a new shelter different from the current shelter.
                newShelter = randi(numShelters);
                while newShelter == currentShelter
                    newShelter = randi(numShelters);
                end

                % Compute entering probability E_p for the new shelter.
                newShelterLeader = shelterLeaders(newShelter);
                newShelterLeaderFitness = fitness_array(newShelterLeader);
                if Fitness > newShelterLeaderFitness
                    E_p = migrationRate * ((Fitness - newShelterLeaderFitness) / Fitness);
                else
                    E_p = 1;  % Candidate is as good as or better than the new shelter leader.
                end

                % Use E_p to decide if the candidate commits to the new shelter.
                if rand < E_p
                    shelters{newShelter} = [shelters{newShelter}, candidate];
 %                   fprintf('         Candidate %d successfully migrated to shelter %d.\n', candidate, newShelter);
                    %pop(candidate,:) = candidateNew;
                    %fitness_array(candidate) = Fitness;
                    break;  % Migration is successful.
                end
                % Otherwise, the loop repeats for another random shelter selection.
            end

            % Step 5. Generalized Search:
            numVars = size(pop, 2);
            % Randomly select two different individuals (ensure they are not candidate itself)
            idx1 = randi(pop_size);
            idx2 = randi(pop_size);
            while idx1 == candidate || idx1 == idx2
                idx1 = randi(pop_size);
            end
            while idx2 == candidate || idx2 == idx1
                idx2 = randi(pop_size);
            end

            % Generate a mutant using a differential-style operation.
            F = 0.5;  % Scaling factor
            mutant = pop(idx1,:) + F*(pop(idx2,:) - pop(idx1,:));
%            fprintf('      Selected individuals %d and %d for mutant generation.\n', idx1, idx2);
%            fprintf('      Mutant vector:\n');
%            disp(mutant);

            % Perform a crossover between the candidate's current solution and the mutant.
            CR = 0.9;  % Crossover probability
            trial = pop(candidate,:);  % Start from the current candidate solution.
            for j = 1:numVars
                if rand < CR
                    trial(j) = mutant(j);
                end
            end
            
%            fprintf('      Candidate %d trial vector after crossover:\n', candidate);
%            disp(trial);

            % Evaluate the new (trial) candidate solution.
            trialFitness = TotalCost(trial, mpc, Nrisk, P_profile, P_renewable);
%            fprintf('      Trial fitness for candidate %d = %.4f (current fitness = %.4f)\n', candidate, trialFitness, fitness_array(candidate));
            if trialFitness < fitness_array(candidate)
%                fprintf('      Mutation accepted for candidate %d. Updating candidate.\n', candidate);
                pop(candidate,:) = trial;
                fitness_array(candidate) = trialFitness;
            end
            % End of explorer branch.


        else
            % Exploiter Branch
 %           fprintf('  -> Candidate %d selected for EXPLOITER branch.\n', candidate);
            % Step 4. Located Search
            if candidate == shelterLeaders(currentShelter)
                currentMembers = shelters{currentShelter};
                % Only perform crossover if there is at least one other candidate in the shelter.
                if length(currentMembers) > 1
                    % Choose a partner candidate from the same shelter (excluding itself)
                    partnerCandidates = currentMembers(currentMembers ~= candidate);
                    partner = partnerCandidates(randi(length(partnerCandidates)));
                    
                    % Perform an arithmetic crossover between the leader and the partner.
                    alpha = rand(); % mixing parameter between 0 and 1
                    candidateCandidate = pop(candidate, :);
                    partnerCandidate   = pop(partner, :);
                    candidateNew = alpha * candidateCandidate + (1 - alpha) * partnerCandidate;
                    
                    % We also add a small mutation to further explore locally.
                    mutationStrength_located = 0.05 * (E_max - E_min);
                    candidateNew = candidateNew + mutationStrength_located * randn(1, numVars);
                    candidateNew = max(candidateNew, E_min);
                    candidateNew = min(candidateNew, E_max);
                    
   %                 fprintf('      Shelter leader candidate %d performing located search with partner %d (alpha = %.2f).\n', candidate, partner, alpha);
   %                 fprintf('         Old candidate vector:\n');
   %                 disp(pop(candidate,:));
   %                 fprintf('         Partner vector:\n');
   %                 disp(partnerCandidate);
   %                 fprintf('         New candidate vector after mutation:\n');
   %                 disp(candidateNew);

                    % Evaluate the new candidate.
                    newFitness = TotalCost(candidateNew, mpc, Nrisk, P_profile, P_renewable);
   %                 fprintf('         New fitness = %.4f (old fitness = %.4f)\n', newFitness, fitness_array(candidate));
                    % Accept the new candidate if it has improved fitness.
                    if newFitness < fitness_array(candidate)
   %                     fprintf('         Located search accepted for candidate %d.\n', candidate);
                        pop(candidate, :) = candidateNew;
                        fitness_array(candidate) = newFitness;
                    end
                end
            else
                % For non-leader exploited candidates, perform a small local mutation.
                mutationStrength_exploiter = 0.05 * (E_max - E_min);
                candidateNew = pop(candidate, :) + mutationStrength_exploiter * randn(1, numVars);
                candidateNew = max(candidateNew, E_min);
                candidateNew = min(candidateNew, E_max);

 %               fprintf('      Expoliting candidate %d performing local mutation.\n', candidate);
 %               fprintf('         Old candidate vector:\n');
 %               disp(pop(candidate,:));
 %               fprintf('         New candidate vector after local mutation:\n');
 %               disp(candidateNew);

                newFitness = TotalCost(candidateNew, mpc, Nrisk, P_profile, P_renewable);
 %               fprintf('         New fitness = %.4f (old fitness = %.4f)\n', newFitness, fitness_array(candidate));
                if newFitness < fitness_array(candidate)
 %                   fprintf('         Local mutation accepted for candidate %d.\n', candidate);
                    pop(candidate, :) = candidateNew;
                    fitness_array(candidate) = newFitness;
                end
            end
        end
    end

    % After processing all candidates, update shelter leaders.
    [shelters, shelterLeaders] = shelterInitialization(fitness_array, numShelters);
    fprintf('\n--- Updated shelter leaders after migration ---\n');
end