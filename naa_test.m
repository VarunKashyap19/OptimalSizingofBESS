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

E_B = NAA_optimizer(Nrisk, mpc, N_BESS, P_profile, P_renewable);