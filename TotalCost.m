function cost = TotalCost(E_B,mpc, Nrisk, P_profile, P_renewable)

[~, ~, success, CO, CD] = compute_charging_schedule(E_B, mpc, Nrisk, P_profile,P_renewable);

if success == 0
    disp(fprintf("Convergence not observed for E_B = %s", num2str(E_B)));
end

CI = computeInvestmentCost(E_B);
CM = computeMaintenanceCost(E_B);

cost = CI + CO + CM + CD;
end

function CM = computeMaintenanceCost(E_B)
    epsilon_b = 20*1000; %per MWh
    CM = epsilon_b * sum(E_B);
end

function CI = computeInvestmentCost(E_B)
    delta_b = 300 *1000; %per MWh
    gamma_b = 24000;
    CI = delta_b * sum(E_B) + gamma_b;
end

