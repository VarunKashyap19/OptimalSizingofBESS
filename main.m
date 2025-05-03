mpc = loadcase('case15da');

n_buses = size(mpc.bus,1);
n_branches = size(mpc.branch,1);

% %Renewable capacities in PU (800 kW = 0.8 MW, 850 kW = 0.85 MW)
solar_capacity = 1; 
wind_capacity_bus11 = 0.65; 
wind_capacity_bus9 = 0.6;

% Update generator data for solar and wind farms
mpc.gen = [
    mpc.gen(:,:); % Existing generator data
    6, 0, 0, 10, -10, 1, 100, 1, solar_capacity, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
    11, 0, 0, 10, -10, 1, 100, 1, wind_capacity_bus11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    9, 0, 0, 10, -10, 1, 100, 1, wind_capacity_bus9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
];

mpc.gencost = [
    mpc.gencost(:,:);
     % New renewable generators
    2   0   0   3   0.01    0.02    0.03;   % Solar farm at Bus 6
    2   0   0   3   0.01    0.03    0.02;   % Wind turbines at Bus 11
    2   0   0   3   0.01    0.025   0.025;  % Wind turbine at Bus 9
];

% Incidence Matrix Generation
E = zeros(n_buses,n_branches);

for k = 1:n_branches
    from_bus = mpc.branch(k,1);
    to_bus = mpc.branch(k,2);

    E(from_bus,k) = 1;
    E(to_bus,k) = -1;
end
E_reduced = E(2:end,:);

%Generating r and x matrix
r_vec = mpc.branch(:,3);
x_vec = mpc.branch(:,4);

r_mat = diag(r_vec);
x_mat = diag(x_vec);


%Creating R and X matrix
R = inv(E_reduced')*r_mat*inv(E_reduced);
X = inv(E_reduced')*x_mat*inv(E_reduced);

%Calculating voltage volatility index
WI = sum(abs(R),2);

%disp("Voltage Volatility Index for each Bus");
%disp(["BUS NUMBER","VOLATGE VOILITION"])
%disp([(2:n_buses)',WI]);


%Algorithm for selecting the appropriate bus
N_BESS = 5;
Nrisk = [];

for i=1:N_BESS
    [~,max_idx] = max(WI);
    selected_bus = max_idx + 1;

    Nrisk = [Nrisk,selected_bus];

    WI = WI - R(:,max_idx);
    WI(max_idx) = 0;
end


%Invocating load profiles
P_profile = zeros((n_buses), 24);
Q_profile = zeros(n_buses, 24);
[P_profile,Q_profile] = load_profile_generator(mpc.bus);

% Extract active and reactive power profiles
P_bus = P_profile(13, :);
Q_bus = Q_profile(13, :);

%Generation profiles
wind_p_pu = [0.65, 0.633, 0.617, 0.600, 0.583, 0.567, 0.550, 0.533, 0.517, 0.500,0.483,0.467,0.450,0.400,0.350,0.300,0.413,0.525,0.638,0.750,0.713,0.675,0.638,0.600];
solar_p_pu = [0,0,0,0,0,0.10,.3 ,.6, .85, .95, 1, 1, 1, 0.95, .8, .6, .3, .1, 0, 0, 0, 0, 0, 0];

wind_q_pu = [0.33, 0.32, 0.31, 0.30, 0.29, 0.28, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.23, 0.20, 0.18, 0.15, 0.21, 0.26, 0.32, 0.38, 0.36, 0.34, 0.32, 0.30];
solar_q_pu = [0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 0.15, 0.30, 0.43, 0.48, 0.50, 0.50, 0.50, 0.48, 0.40, 0.30, 0.15, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00];


wind_bus9 = wind_capacity_bus9 * wind_p_pu;
wind_bus11 = wind_capacity_bus11 *  wind_p_pu;
solar_bus6 = solar_capacity * solar_p_pu;


%Making renewable power matrix
P_renewable = zeros(n_buses, 24);
P_renewable(6, :) = solar_bus6;
P_renewable(9, :) = wind_bus9;
P_renewable(11, :) = wind_bus11;

wind_q_bus9 = wind_capacity_bus9 * wind_q_pu;
wind_q_bus11 = wind_capacity_bus11 * wind_q_pu;
solar_q_bus6 = solar_capacity * solar_q_pu;

% Create a Q renewable generation matrix (similar to P_renewable)
Q_renewable = zeros(n_buses, 24);
Q_renewable(6, :) = solar_q_bus6;
Q_renewable(9, :) = wind_q_bus9;
Q_renewable(11, :) = wind_q_bus11;



%Creating Voltage profiles
time_horizon = 24;

voltage_profiles = zeros(n_buses-1, time_horizon);
V_slack = 1.02;  % Slack bus voltage in p.u.

for h = 1:time_horizon
    P_net = zeros(n_buses-1,1);
    Q_net = zeros(n_buses-1,1);
    for i = 2:n_buses
        idx = i-1;
        load_p = P_profile(i,h);
        load_q = Q_profile(i,h);
        
        gen_p = 0;
        if i == 6
            gen_p = solar_bus6(h);  % Solar generation at bus 6
        elseif i == 9
            gen_p = wind_bus9(h);   % Wind generation at bus 9
        elseif i == 11
            gen_p = wind_bus11(h);  % Wind generation at bus 11
        end
        
        P_net(idx) = gen_p - load_p;
        Q_net(idx) = -load_q;
    end
    % Compute voltage deviations for buses 2-15
    dV = R * P_net + X * Q_net; 
    % Bus voltage = slack voltage plus deviation.
    voltage_profiles(:,h) = V_slack + dV;
end

voltage_profiles = [1.02 * ones(1, 24); voltage_profiles];

% Plotting Voltage Profiles using DistFlow
figure;
hold on;
for i = 1:(time_horizon)
    plot(1:n_buses, voltage_profiles(:,i), 'DisplayName', ['Hour =  ' num2str(i)]);
end
xlabel('Buses');
ylabel('Voltage (p.u.)');
title('Voltage Profiles for IEEE 15-Bus System (Buses 2-15)');
ylim([.90,1.05]);
grid on;
hold off;

%Comparing with AC Load flow
power_flow_results = cell(1, time_horizon);
mpc.bus(1,8) = 1.02;

for h = 1:time_horizon
    mpc_h = mpc;

    mpc_h.gen(1,6) = 1.02;

    % Update loads
    for i = 1:n_buses
        mpc_h.bus(i, 3) = P_profile(i, h);
        mpc_h.bus(i, 4) = Q_profile(i, h);
    end
    
    % Update renewable generation
    mpc_h.gen(end-2, 2) = solar_bus6(h);
    mpc_h.gen(end-1, 2) = wind_bus11(h);
    mpc_h.gen(end, 2) = wind_bus9(h);
    
    results = runpf(mpc_h, mpoption('PF_ALG', 1, 'VERBOSE', 0,'OUT_ALL', 0));

    voltage_profiles(:, h) = results.bus(:, 8);
    power_flow_results{h} = results;
end

% Display and plot results
figure;
plot(1:n_buses, voltage_profiles');
yline(0.97, 'r--', 'LineWidth', 1.5); 
yline(1.03, 'r--', 'LineWidth', 1.5);
xlabel('Buses');
ylabel('Voltage (p.u.)');
ylim([.90,1.05]);
title('24-Hour AC Voltage Profiles for IEEE 15-Bus System');
legend(cellstr(num2str((1:n_buses)')), 'Location', 'eastoutside');
% Count how many violations occurred
violations = sum(sum(voltage_profiles < 0.97 | voltage_profiles > 1.03));
disp(['Number of voltage violations: ', num2str(violations), ' out of ', num2str(n_buses * time_horizon), ' total measurements']);
grid on;




% Calling NAA Optimizer
E_B = NAA_optimizer(Nrisk, mpc, N_BESS, P_profile, P_renewable);
% Display optimization results
disp('Final BESS Capacities:');
disp(E_B);



% Generating p and q achedule for the optimized BESS capacities
[p_schedule, q_schedule, success, CO, CD] = compute_charging_schedule(E_B, mpc, Nrisk, P_profile,P_renewable);



%Generation profiles plot (can make both wind profiles different)
figure;
plot(wind_bus9, 'DisplayName', 'Wind Active Power Generation at bus 9');
hold on;
plot(wind_bus11, 'DisplayName', 'Wind Active Power Generation at bus 11');
plot(solar_bus6, 'DisplayName', 'Solar Active Power Generation at bus 6');
plot(wind_q_bus9, 'DisplayName', 'Wind Reactive Power Generation at bus 9');
plot(wind_q_bus11, 'DisplayName', 'Wind Rective Power Generation at bus 11');
plot(solar_q_bus6, 'DisplayName', 'Solar Rective Power Generation at bus 6');
hold off;
ylim([0, 1.2]);
legend('show');
xlabel('Hour');
ylabel('Generation in MW/MVar');
title('Generation profiles');

% Plot active and reactive power profiles
figure;
plot(1:24, P_bus, 'DisplayName', 'Active Power (P)');
hold on;
plot(1:24, Q_bus, 'DisplayName', 'Reactive Power (Q)');
hold off;

xlabel('Hour');
ylabel('Power (MW/MVar)');
title('Active and Reactive Power Profile for Bus 13');
legend('Location', 'best');
grid on;




% Create a figure for BESS charging/discharging schedule
figure;

% Setup time vector (24 hours)
time = 1:24;

% Plot BESS active power schedule (p_schedule)
hold on;
colors = lines(length(Nrisk));
legendEntries = cell(length(Nrisk) + 1, 1);

for i = 1:length(Nrisk)
    plot(time, p_schedule(i,:), 'LineWidth', 2, 'Color', colors(i,:));
    legendEntries{i} = sprintf('BESS at Bus %d', Nrisk(i));
end

% Plot total power from all BESS units
total_power = sum(p_schedule, 1);
plot(time, total_power, 'k--', 'LineWidth', 2);
legendEntries{end} = 'Total BESS Power';

% Add horizontal line at y=0 to differentiate charging from discharging
yline(0, 'k-', 'LineWidth', 1);

% Add labels and formatting
grid on; grid minor;
xlabel('Time (hour)', 'FontSize', 12);
ylabel('Active Power (MW)', 'FontSize', 12);
title('BESS Charging/Discharging Schedule', 'FontSize', 14);
legend(legendEntries, 'Location', 'best', 'FontSize', 10);

% Add charging/discharging labels
max_power = max(abs(p_schedule(:)));
text(2, max_power*0.7, 'Discharging (Power > 0)', 'FontWeight', 'bold');
text(2, -max_power*0.7, 'Charging (Power < 0)', 'FontWeight', 'bold');

% Set y-axis limits to be symmetric around zero
ylim([-max_power*1.1, max_power*1.1]);

% Add a second figure for State of Charge
figure;
hold on;

% Get SOC values that were created in compute_charging_schedule
soc = evalin('base', 'soc');

% Plot SOC for each BESS
for i = 1:length(Nrisk)
    plot(time, soc(i,:), 'LineWidth', 2, 'Color', colors(i,:));
end

% Add formatting
grid on; grid minor;
xlabel('Time (hour)', 'FontSize', 12);
ylabel('State of Charge', 'FontSize', 12);
title('BESS State of Charge', 'FontSize', 14);
legend(legendEntries(1:end-1), 'Location', 'best', 'FontSize', 10);
ylim([0, 1]);

% Add horizontal lines for SOC limits
SOC_min = 0.2;
SOC_max = 0.9;
yline(SOC_min, 'r--', 'LineWidth', 1);
yline(SOC_max, 'r--', 'LineWidth', 1);
text(1, SOC_min-0.05, 'SOC_{min}', 'FontWeight', 'bold');
text(1, SOC_max+0.05, 'SOC_{max}', 'FontWeight', 'bold');



p_inj = sum(P_profile,1) - sum(P_renewable,1);
p_demand = zeros(15,24);

for t = 1:24    
    for i = 2:n_buses
                        [isInNrisk, indexInNrisk] = ismember(i, Nrisk); % Check membership and get index

                        if isInNrisk
                            p_demand(i,t) = (P_profile(i,t) - P_renewable(i,t) - p_schedule(indexInNrisk,t));
                        else
                            p_demand(i,t) = (P_profile(i,t) - P_renewable(i,t));
                        end
     end
end

p_demand = sum(p_demand,1);
figure
plot(1:24,p_inj,'DisplayName', 'Injected power without bess');
hold on
plot(1:24,p_demand,'DisplayName', 'Injected power with optimal capacity of bess');
legend('show');

voltage_profiles_without_bess = voltage_profiles;
voltage_profiles_with_bess = zeros(n_buses, time_horizon);

for h = 1:time_horizon
    mpc_h = mpc;
    mpc_h.gen(1,6) = 1.02;  % Slack bus voltage

    % Update loads
    for i = 1:n_buses
        mpc_h.bus(i, 3) = P_profile(i, h);  % Active power load
        mpc_h.bus(i, 4) = Q_profile(i, h);
    end
    
    % Update renewable generation
    mpc_h.gen(end-2, 2) = solar_bus6(h);    % Solar at bus 6
    mpc_h.gen(end-1, 2) = wind_bus11(h);    % Wind at bus 11
    mpc_h.gen(end, 2) = wind_bus9(h);       % Wind at bus 9
    
    % Add BESS power injections to the buses where BESS is installed
    for i = 1:length(Nrisk)
        bess_bus = Nrisk(i);
        mpc_h.bus(bess_bus, 3) = mpc_h.bus(bess_bus, 3) - p_schedule(i, h);  
        mpc_h.bus(bess_bus, 4) = mpc_h.bus(bess_bus, 4) - q_schedule(i, h);  
    end
    
    % Run power flow with BESS
    results = runpf(mpc_h, mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0));
    
    voltage_profiles_with_bess(:, h) = results.bus(:, 8);
end

% Display and plot results with BESS
disp('Voltage Profiles with BESS (p.u.) for all buses over 24 hours:');
disp(voltage_profiles_with_bess);

% Plot comparing voltage profiles before and after BESS
figure;
plot(1:n_buses, voltage_profiles_with_bess');
yline(0.97, 'r--', 'LineWidth', 1.5); 
yline(1.03, 'r--', 'LineWidth', 1.5);
xlabel('Buses');
ylabel('Voltage (p.u.)');
ylim([0.90, 1.05]);
title('24-Hour AC Voltage Profiles WITH optimal capacity of BESS');
grid on;

% Check if voltage profiles are within limits
voltage_min = min(voltage_profiles_with_bess(:));
voltage_max = max(voltage_profiles_with_bess(:));

disp(['Minimum voltage with BESS: ', num2str(voltage_min), ' p.u.']);
disp(['Maximum voltage with BESS: ', num2str(voltage_max), ' p.u.']);

if voltage_min >= 0.97 && voltage_max <= 1.03
    disp('Voltage profile is within acceptable limits (0.97-1.02 p.u.)');
else
    disp('Voltage profile exceeds acceptable limits!');
    
% Count how many violations occurred
violations = sum(sum(voltage_profiles_with_bess < 0.97 | voltage_profiles_with_bess > 1.03));
disp(['Number of voltage violations: ', num2str(violations), ' out of ', num2str(n_buses * time_horizon), ' total measurements']);
end