function [p_schedule, q_schedule,success, CO, CD] = compute_charging_schedule_with_q(E_B, mpc, Nrisk, P_profile, P_renewable, Q_profile, Q_renewable, R, X)

    E_B = E_B(:);

    % Dimensions
    N_buses = size(mpc.bus, 1); 
    N_time = 24;  
    N_bess = length(E_B);
    
    % Parameters (example values)
    days = 365;
    eta_c = 0.05;
    delta_b = 300*1000; %per MWh 
    gamma_b = 24000;
    d = 0.05;
    SOC_min = 0.2;
    SOC_max = 0.9;
    Tb = 1000;
    kB = 0.95;
    a_gen = 1875;
    b_gen = 0;
    c_gen = 0;
    V_min = 0.97;
    V_max = 1.03;

    % Initialize outputs
    p_schedule = zeros(N_bess, N_time);
    q_schedule = zeros(N_bess, N_time);
    success = false;
    
    
    % Create the bus admittance matrix (in per-unit)
    [Ybus, ~, ~] = makeYbus(1, mpc.bus, mpc.branch);
    

    try
        cvx_begin quiet
            variables p(N_bess, N_time) q(N_bess, N_time) soc(N_bess, N_time)
            variables W(N_buses, N_buses, N_time)

            expressions P_injection(N_buses, N_time) Q_injection(N_buses, N_time) Ybus_m(N_buses, N_buses) P_profile_m(N_buses, N_time) P_renewable_m(N_buses, N_time) pbess_total(N_bess, 1)
            expressions R_m((N_buses-1), (N_buses-1)) X_m((N_buses-1), (N_buses-1)) V_slack Q_profile_m(N_buses, N_time) Q_renewable_m(N_buses, N_time) qbess_total(N_bess, 1)
            
            Ybus_m = full(Ybus);
            P_profile_m = P_profile;
            P_renewable_m = P_renewable;
            Q_profile_m = Q_profile;
            Q_renewable_m = Q_renewable;

            %expression P1_temp(24)
            for t = 1:N_time
                W_t = W(:,:,t); 
                for i = 1:N_buses
                    P_injection(i,t) = real(sum(conj(Ybus_m(i,:)) .* W_t(i,:)));
                    Q_injection(i,t) = imag(trace(conj(Ybus_m(i,:)) * W_t(:,i)));
                end
            end

            %Calculating Voltage Violitions
            expressions voltage_profiles((N_buses-1), 24) violitions
            V_slack = 1.02;
            R_m = R;
            X_m = X;
            dV = R_m * (-1 * P_injection(2:end,:)) + X_m *(-1 * Q_injection(2:end,:)); 
            % Bus voltage = slack voltage plus deviation.
            voltage_profiles = V_slack + dV;
            violations = 0;

            for i = 1:(N_buses-1)
                for t = 1:N_time
                    violations = violations + max(0, 0.97 - voltage_profiles(i,t));
                    violations = violations + max(0, voltage_profiles(i,t) - 1.03);
                end
            end

            expression P1(1,24)
            P1 = P_injection(1, :);
            cost_gen = a_gen * (P1.^2) + b_gen * P1 + c_gen;
            cost_BESS = (1 - eta_c) * sum(abs(p), 1);  % returns a 1xN_time vector
            expression operational_cost
            operational_cost = sum(cost_gen + cost_BESS) * days;
            LC = 2 * E_B * Tb;  
            LC_rest = LC - days * sum(abs(p), 2);  
            expression residual_cost
            residual_cost = - sum( ((delta_b) .* E_B .* (LC_rest ./ LC)) / (1 + d) + gamma_b .* kB );
            expression violition_cost
            v_c = 2e8;
            violation_cost = violations*v_c;

            % Objective: minimize total cost.
            minimize( operational_cost + residual_cost + violation_cost);

            %Defining Constraints

            subject to

                %Slack Bus Constraint
                W(1,1,:) == 1.02;

                %Cone Constraints
                for t = 1:N_time
                    for k = 1:size(mpc.branch, 1)
                        i = mpc.branch(k, 1);
                        j = mpc.branch(k, 2);                       
                        norm([2*W(i,j,t); W(i,i,t)-W(j,j,t)], 2) <= W(i,i,t) + W(j,j,t);
                       
                    end
                end

                %SOC Constraints
                SOC_init = ones(N_bess, 1) * 0.5; % Example: start at 50% SOC
                for i = 1:N_bess
                    soc(i,1) == SOC_init(i);
                    for t = 1:(N_time-1)
                        % SOC dynamics equation (assuming 1-hour time steps)
                        soc(i,t+1) == soc(i,t) - p(i,t)/E_B(i);
                    end
                    % SOC bounds for all time periods
                    soc(i,:) >= SOC_min;
                    soc(i,:) <= SOC_max;
                end



                %Power Balance Constraint
                expression P_demand(1,24)
                P_demand = sum(P_profile_m,1) - sum(P_renewable_m,1) - sum(p,1);
                
                for i = 1:N_time
                    P_injection(1,i) >= P_demand(1,i)
                end


                for t = 1:N_time
                    for i = 2:N_buses
                        [isInNrisk, indexInNrisk] = ismember(i, Nrisk); % Check membership and get index

                        if isInNrisk
                            P_injection(i,t) == (P_profile_m(i,t) - P_renewable_m(i,t) - p(indexInNrisk,t));
                        else
                            P_injection(i,t) == (P_profile_m(i,t) - P_renewable_m(i,t));
                        end
                    end
                end
                expression Q_demand(1,24)
                Q_demand = sum(Q_profile_m,1) - sum(Q_renewable_m,1) - sum(q,1);
                
                for i = 1:N_time
                    Q_injection(1,i) >= Q_demand(1,i)
                end


                for t = 1:N_time
                    for i = 2:N_buses
                        [isInNrisk, indexInNrisk] = ismember(i, Nrisk); % Check membership and get index

                        if isInNrisk
                            Q_injection(i,t) == (Q_profile_m(i,t) - Q_renewable_m(i,t) - q(indexInNrisk,t));
                        else
                            Q_injection(i,t) == (Q_profile_m(i,t) - Q_renewable_m(i,t));
                        end
                    end
                end


                %Voltage Constraint
                for t = 1:N_time
                    for i = 1:N_buses
                        W(i,i,t) >= (V_min)^2;
                        W(i,i,t) <= (V_max)^2;
                    end
                end

                %Power balance constraint
                for i = 1:N_bess
                    sum(p(i,:)) == 0; % Net zero energy constraint for each BESS
                end
                
                %Battery Operating Constraints
                P_B = E_B / 2;
                for i = 1:N_bess
                    for t = 1:N_time
                        norm([p(i,t); q(i,t)], 2) <= P_B(i);
                    end
                end

                

        cvx_end

        if strcmp(cvx_status, 'Solved')
            success = true;
            p_schedule = p;
            q_schedule = q;
            CO = operational_cost;
            CD = residual_cost;
        else
            warning('Optimization did not converge.');
        end

    catch ME
        disp('Error during optimization:');
        disp(ME.message);
    end

    assignin("base","W", W);
    assignin("base","p", P_injection);
    assignin("base","Ybus", Ybus_m);
    assignin("base","soc", soc);
    assignin("base",'violations_fxn',violation_cost);

end