%%%%%
% MATLAB function to run main solver loop for UNNAMED ECON IAM SOLVER THINGY. 
% NOTES: 
%   - The equilibrium controller uses parallelization for lsqnonlin calls.
%   - Applying the launch mask to more locations generally but not always makes the equilibrium solver run faster.
%   - Parallel pool is recreated and closed each time the equilibrium controller is called, in openAccessSolver. Probably best to pull this out and put it in the main script, but it is nice to be able to call openAccessSolver with parallelization independently in one line.
%   - NOTES ON PARALLELIZATION:
%       - The equilibrium controller uses parallelization for lsqnonlin calls, parallelizing over altitude shells. Default setting is to use 40 cores since there are 40 altitude shells. This parameter is set from the bash script conductor.sh as a user-configurable input.
%       - The GMPHD propagator assumes debris are continuous over altitude and radius. For reporting, it integrates over the size and altitude distributions to calcualte numbers of debris of three types in all altitude shells. The post-processing function PostprocessDebris uses parallelization for Gaussian integration over time steps. Default setting is to use 10 cores since there are 10 time steps. This parameter is NOT set from the bash script conductor.sh as a user-configurable input -- ctrl-F for PostprocessDebris and adjust last argument to change this.
%
% Author: Akhil Rao
% Last updated: September 2023
%%%%%
function iam_solver(stem, model_type, launch_pattern_type, parameter_file, model_horizon_str, n_workers)

big_timer = tic; % Start the timer

%%%% SIMULATION/SCENARIO PARAMETERS %%%%
% %%% BEGIN block of inputs for testing
% stem = 'destroying-termite-verbalize-Zinnia'; %Stem for output files
% model_type = 'GMPHD'; %Type of physical model to use
% launch_pattern_type = 'equilibrium'; %Type of launch pattern to use
% parameter_file = 'scenarios/parsets/5yr-rule-0.5-Compliance-parset.csv'; %Parameter file to use. Use 'benchmark' for default parameters
% model_horizon_str = '10'; %Length of simulation in years
% n_workers = '40'; % Number of workers for parallelization in equilibrium case. Set to 1 for serial. Check this for your system.
%%% END testing block

model_horizon = str2double(model_horizon_str); % Convert char input to double
tf = [0:model_horizon]; %Years of simulation
filename = strcat(stem, '-', model_type, '-', launch_pattern_type, '-', model_horizon_str,'yrs__out-data', '.csv'); %Output filename
n_workers = str2double(n_workers); % Convert char input to double

% Print stem
fprintf('Stem: %s\n', stem);

%% Add path to MOCAT4S and GMPHD filter files
addpath('MOCAT4S')
addpath('GMPHD')

%%%% CONSTELLATION BLOCK
% Set up constellation parameters. Read in from constellation-parameters.csv in parsets folder. Change this to run different constellation scenarios.
constellation_params = constellation_parameters('scenarios/parsets/constellation-parameters.csv');
%%%% END CONSTELLATION BLOCK

%%%% LOCATION MASKING BLOCK
% TODO: Make this a switch based on string input. Definition also appears in scenarioNamer.m
% Create mask(s): 1 if location is under open access, 0 if access is restricted/blocked. Naming convention is [where it's 1]-mask. Useful to study keep-out zones.
launch_mask = ones(1,40);
%%%% END LOCATION MASKING BLOCK

% Set up physical model parameters from MOCAT
VAR = MOCAT4S_VAR_Cons();

% Set up physical model parameters from GMPHD
GMPHD_params = GMPHD_VAR_Cons();

% Set up economic model parameters
econ_params = set_econ_parameters(VAR);

%%%% PHYSICAL MODEL-SPECIFIC CODE BLOCKS %%%%
%%%% BEGIN MOCAT-4S LOADING BLOCK %%%%
if strcmp(model_type, 'MOCAT')
    % Read in initial populations from TLE data 
    filename_N = 'x0_TLE/Counts_DEBRIS_bins_35.csv';
    filename_S = 'x0_TLE/Counts_PAYLOADslot_bins_35.csv';
    filename_Su = 'x0_TLE/Counts_PAYLOADunslot_bins_35.csv';
    filename_D = 'x0_TLE/Counts_DERELICT_bins_35.csv';
    Ni = readmatrix(filename_N);
    Ni = Ni(:,2)';  % Extract initial condition for debris
    Si = readmatrix(filename_S);
    Si = Si(:,2)';  % Extract initial condition for slotted payloads
    Sui = readmatrix(filename_Su);
    Sui = Sui(:,2)';  % Extract initial condition for unslotted payloads
    Di = readmatrix(filename_D);
    Di = Di(:,2)'; % Extract initial condition for derelict objects
    % Set initial conditions
    x0 = [Si';Di';Ni';Sui'];
end
%%%% END MOCAT-4S LOADING BLOCK %%%%

%%%% BEGIN GMPHD LOADING BLOCK %%%%
%%%%%% NOTE: GMPHD model currently does not load initial conditions from TLE data. This is a TODO for future development.
if strcmp(model_type, 'GMPHD')
    ntime = 10;  % Number of integration time steps within year
    nyears = 1;
    dicoll=10; % Number of timesteps between collisions
    % Set up initial conditions
    %% NULL INITIAL CONDITION BLOCK -- ALL ZEROS. Useful for testing.
    % constellation_levels = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    %                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    %                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    %                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    % fringe_levels = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    %                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    %                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    %                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    %% TLE INITIAL CONDITION BLOCK
    filename_S = 'x0_TLE/Counts_PAYLOADslot_bins_35.csv';
    filename_Su = 'x0_TLE/Counts_PAYLOADunslot_bins_35.csv';
    constellation_levels = readmatrix(filename_S);
    constellation_levels = constellation_levels(:,2)';  % Extract initial condition for slotted payloads
    fringe_levels = readmatrix(filename_Su);
    fringe_levels = fringe_levels(:,2)';  % Extract initial condition for unslotted payloads  
    % Initialize launches as replacement rate on initial stocks
    constellation_launches = 1/(+VAR.Dt) * constellation_levels;
    fringe_launches = 1/(+VAR.Dt) * fringe_levels;
    % Set up GMPHD parameters
    GMPHD_params.r0 = 700+GMPHD_params.rE; % reference height for drag
    % Build shellstruct to hold discretized objects for reporting and collision calculations
    shellstruct = buildShellStruct(constellation_levels, fringe_levels, constellation_launches, fringe_launches, GMPHD_params);
    % Build initial GM objects
    initial_GM_objects = prepareInitialGMs(GMPHD_params, nyears, ntime, dicoll);
    t = initial_GM_objects.t;
    GM = initial_GM_objects.GM;
    i_collisions = initial_GM_objects.i_collisions;
    GM_coll = initial_GM_objects.GM_coll;
    % Create initial debris counts. This involves propagating the model one time step with zero launch rates, then using the PostprocessDebris function to calculate sizes of three debris types: LNTs (0-10cm), small trackable fragments (10-50cm), and large trackable fragments (50cm+).
    % NOTE: this section parallelizes Gaussian integration over the time steps used for propagator integration. So if t = 10, the code assumes you have 10 cores available for parallelization. CHECK THIS FOR YOUR SYSTEM.
    %% Propagation step to generate GM1
    tic
    [t,GM1,shellstruct] = PropagateGM_timehist_collisions(t, GM, GMPHD_params, shellstruct, i_collisions, GM_coll);
    toc
    %% Define bounds of integration for Gaussians. This is used throughout the rest of the code.
    rbounds=[0, 0.1/1000, 0.5/1000, 3/1000]; % Size categories for debris in post-processing. 0-10cm, 10-50cm, 50cm-1m
    abounds=[200:35:1600]+GMPHD_params.rE; % altitude bins for debris in post-processing. 200-1600 km in 35 km blocks to match MOCAT for prototype.
    %% Numerical integration over Gaussians. Runs in parallel.
    fprintf('Now starting debris post-processing...\n');
    tic
    [D_lnt, D_st, D_lt, Si, Sui, lam]=PostprocessDebris(GM1, t, abounds, rbounds, GMPHD_params, shellstruct, 0, min(length(t),40) );
    toc
    x0 = [Si';D_lt';D_st';D_lnt';Sui'];
end
%%%% END GMPHD LOADING BLOCK %%%%

%% Write empty file to folder if this is the benchmark case
if strcmp(parameter_file, 'benchmark')
    fid = fopen(strcat('scenarios/', stem, '/', parameter_file), 'w');
    fclose(fid);
end

%%%% UPDATE PARAMETERS FOR SCENARIOS IF SUPPLIED
% Modify parameters if a specific parameter file is provided (not "benchmark"), then write parameter file to output folder
if ~strcmp(parameter_file, 'benchmark')
    [VAR, econ_params, constellation_params] = modifyParameters(VAR, econ_params, constellation_params, parameter_file);
    fprintf('Parameters modified\n');
    % Strip "scenarios/parsets/" from string in "parameter_file", and write an empty file with that name into the output folder
    parameter_file_stem = extractAfter(parameter_file, "scenarios/parsets/");
    parameter_file_stem = extractBefore(parameter_file_stem, ".csv");
    fid = fopen(strcat('scenarios/', stem, '/', parameter_file_stem), 'w');
    fclose(fid);
end
%%%% END PARAMETER UPDATE BLOCK

%% Build the cost function
cost_fn_params = buildCostFunction(VAR, econ_params);

fprintf('Cost function parameters:\n', 0); % Print the cost function parameters
cost_fn_params.k_star
cost_fn_params.cost
%% Testing -- check components in detail
% cost_fn_params.total_lift_price;
% cost_fn_params.deorbit_maneuver_cost';
% cost_fn_params.stationkeeping_cost';
% cost_fn_params.lifetime_loss_cost';
% cost_fn_params.v_drag;

%% Push k_star to VAR after any disosalTime modifications have been made, push k_star and cost parameters into econ_params
VAR.k_star = cost_fn_params.k_star;
econ_params.k_star = cost_fn_params.k_star;
econ_params.total_lift_price = cost_fn_params.total_lift_price;
econ_params.deorbit_maneuver_cost = cost_fn_params.deorbit_maneuver_cost;
econ_params.stationkeeping_cost = cost_fn_params.stationkeeping_cost;
econ_params.lifetime_loss_cost = cost_fn_params.lifetime_loss_cost;
econ_params.cost = cost_fn_params.cost;

%% Integrate over multiple time periods with MOCAT-4S default feedback controller yearly determined launch rate
% Assign initial condition to X for accumulating states
X = x0';

% Create location indices
location_indices = linspace(1,40,40);

%%%% INITIAL PERIOD LAUNCH RATE
lam(:,1) = 1/(+VAR.Dt) * Si'; % Slotted objects -- approximate replacement rate feedback rule
% CONSTELLATION BUILDER
%% Call the constellationBuildup function to get the modified launch rate given the constellation buildup parameters. Once for each constellation.
for i = 1:constellation_params.n_constellations
    lam(constellation_params.location_index(i),1) = constellationBuildup(constellation_params.location_index(i), constellation_params.final_size(i), constellation_params.linear_rate(i), Si');
end

% Generate initial period launch rate
tic; % Start the timer
fprintf('Now starting period %d...\n', 0); % Print the starting message

%% Fringe population autonomous controller. 
% TODO---profiles need to be created and tested. Code block describes intent for future development.
% if strcmp(launch_pattern_type, 'autonomous')
%     % Using a launch rate profile over shells that is read in from a file. The file should be a CSV with a header row and 1 + T columns: 1 shell number and T launch rates for T periods (0, 1, ..., T).
%     % Read in the launch rate profile
%     launch_rate_profile = readmatrix('launch-rate-profiles/launch-rate-profile-425-1250.csv');
%     lam(:,2) = launch_rate_profile(:,2); % Unslotted objects launched in period 0
% end

%% Fringe sat-pop feedback controller
if strcmp(launch_pattern_type, 'sat_feedback')
    % Using 5% of shell satellite population as launch rate
    % lam(:,2) = 0.05 * Sui'; % Unslotted objects -- 5% feedback rule
    lam(:,2) = 1/(+VAR.Dt) * Sui' .* launch_mask'; % Unslotted objects -- approximate replacement rate feedback rule
end

%% Fringe equilibrium controller
if strcmp(launch_pattern_type, 'equilibrium')
% Using open-access equilibrium launch rate

    % Create guess
    solver_guess = [0.05 * Sui'.*launch_mask']';
    lam(:,2) = solver_guess; % Unslotted objects -- 5% feedback rule
    
    % MOCAT block
    if strcmp(model_type, 'MOCAT')
        %% Generate MOCAT output
        tspan = linspace(0,10,VAR.N_step);
        OUT = MOCAT4S(tspan, x0, lam, VAR);
        %% Testing -- check econ values
        % ror = fringeRateOfReturn("linear", econ_params, OUT, location_indices, launch_mask);
        % for i = 1:VAR.N_shell
        %     collision_probability(i) = calculateCollisionProbability_MOCAT4S(OUT, VAR, i);
        % end
        % Solve for equilibrium launch rates
        lam(:,2) = openAccessSolver(solver_guess, launch_mask, "MOCAT", VAR, x0, "linear", econ_params, lam(:,1), n_workers)'; % Unslotted objects = Fringe satellites
    end

    % GMPHD block
    if strcmp(model_type, 'GMPHD')
        % Push lam back into shellstruct launch rates
        for i = 1:40 % Hard-coded to match MOCAT shell discretization for launch and reporting.
            shellstruct(i).launch_rate_const(end) = lam(i,1) ;
            shellstruct(i).launch_rate_fringe(end) = lam(i,2);
        end
        % Propagate
        GM1_last = GM1(end,:); %% Extract last element of GM1 to use as initial condition for next period
        [t,GM1,shellstruct] = PropagateGM_timehist_collisions(t, GM1_last, GMPHD_params, shellstruct, i_collisions, GM_coll);
        % Extract for reporting
        [D_lnt, D_st, D_lt, Si, Sui, lam]=PostprocessDebris(GM1, t, abounds, rbounds, GMPHD_params, shellstruct, 0, length(t));
        % Package all the GMPHD parameters into a single struct to pass to openAccessSolver
        GMPHD_all.t = t;
        GMPHD_all.GMPHD_params = GMPHD_params;
        GMPHD_all.shellstruct = shellstruct;
        GMPHD_all.GM_coll = GM_coll;
        GMPHD_all.i_collisions = i_collisions;
        % Solve for equilibrium launch rates
        lam(:,2) = openAccessSolver(solver_guess, launch_mask, 'GMPHD', GMPHD_all, GM1, "linear", econ_params, lam(:,1), n_workers)';
    end
end
elapsedTime = toc; % Stop the timer and get the elapsed time
fprintf('Time taken for period 0: %.2f seconds\n', elapsedTime); % Print the time taken

% Assign initial condition to LAM for accumulating launch rates
LAM1 = lam(:,1)';
LAM2 = lam(:,2)';

%%% BEGIN LITTLE TEST CODE
%% Useful for testing openAccessSolver
% tic;
% init_launch_rate = openAccessSolver(solver_guess, launch_mask, "MOCAT", VAR, x0, "linear", econ_params, lam(:,1), 40);
% elapsedTime = toc; % Stop the timer and get the elapsed time
% fprintf('Time taken for test: %.2f seconds\n', elapsedTime); % Print the time taken
%%% END LITTLE TEST CODE

%%%% REST OF THE SIMULATION

% Main loop over time periods
for timeIDX = 1:length(tf)-1
    tic; % Start the timer
    fprintf('Now starting period %d...\n', timeIDX); % Print the starting message

    %%%% BEGIN MOCAT-4S PROPAGATOR
    if strcmp(model_type, 'MOCAT')
        tspan = linspace(tf(timeIDX),tf(timeIDX+1),VAR.N_step); % Create linearly-spaced vector from tf(t):tf(t+1) with VAR.N_step points. This is the within-year time-stepper for integration.
        OUT = MOCAT4S(tspan, x0, lam, VAR); % Propagate for the first time period
        lam(:,1) = 1/(+VAR.Dt) * OUT.S(end,:)';  % approximate replacement rate feedback rule
    end
    %%%% END MOCAT-4S PROPAGATOR

    %%%% BEGIN GMPHD PROPAGATOR
    if strcmp(model_type, 'GMPHD')
        % Push lam back into shellstruct launch rates
        for i = 1:40 % Hard-coded to match MOCAT shell discretization for launch and reporting.
            shellstruct(i).launch_rate_const(end) = lam(i,1) ;
            shellstruct(i).launch_rate_fringe(end) = lam(i,2);
        end
        % Propagate
        GM1_last = GM1(end,:); %% Extract last element of GM1 to use as initial condition for next period
        [t,GM1,shellstruct] = PropagateGM_timehist_collisions(t, GM1_last, GMPHD_params, shellstruct, i_collisions, GM_coll);
        % Extract for reporting
        [D_lnt, D_st, D_lt, Si, Sui, lam]=PostprocessDebris(GM1, t, abounds, rbounds, GMPHD_params, shellstruct, 0, length(t));
        % Package all the GMPHD parameters into a single struct to pass to openAccessSolver
        lam(:,1) = 1/(+VAR.Dt) * Si';  % approximate replacement rate feedback rule
    end
    %%%% END GMPHD PROPAGATOR

    % Launch rates for next time step based on states in previous time step
    %% Constellation/slotted satellite controller. 
    % lam(:,1) = 0.05 * OUT.S(end,:)';  % 5% feedback rule
    % Call the constellationBuildup function to get the modified launch rate given the constellation buildup parameters. Once for each constellation.
    for i = 1:constellation_params.n_constellations
        lam(constellation_params.location_index(i),1) = constellationBuildup(constellation_params.location_index(i), constellation_params.final_size(i), constellation_params.linear_rate(i), Si');
    end

    %% Fringe population feedback controller
    if strcmp(launch_pattern_type, 'sat_feedback')
    % Using 5% of shell satellite population as launch rate
        tspan = linspace(tf(timeIDX),tf(timeIDX+1),VAR.N_step); % Create linearly-spaced vector from tf(t):tf(t+1) with VAR.N_step points. This is the within-year time-stepper for integration.
        OUT = MOCAT4S(tspan, x0, lam, VAR); % Propagate for the first time period
        lam(:,2) = 1/(+VAR.Dt) * OUT.Su(end,:)'; % Unslotted objects, approximate replacement rate feedback rule
    end

    %% Fringe equilibrium controller
    if strcmp(launch_pattern_type, 'equilibrium')
    % Using open-access equilibrium launch rate

        %% MOCAT block
        if strcmp(model_type, 'MOCAT')

            % Create solver guess
            ror = fringeRateOfReturn("linear", econ_params, OUT, location_indices, launch_mask, 'MOCAT');
            for i = 1:VAR.N_shell
                collision_probability(i) = calculateCollisionProbability_MOCAT4S(OUT, VAR, i);
            end

            %% Generate econ solver guess
            solver_guess = lam(:,2)' - lam(:,2)'.*(ror - collision_probability).*launch_mask; % Previous period's launch rate works well near steady state. Add excess return scaled launch rate to perturb in the direction of excess return.
            fprintf('Solver guess created\n');
            solver_guess

            % Solve for equilibrium launch rates
            lam(:,2) = openAccessSolver(solver_guess, launch_mask, "MOCAT", VAR, x0, "linear", econ_params, lam(:,1), n_workers)'; % Unslotted objects = Fringe satellites
            fprintf('Equilibrium condition solved\n');
            lam(:,2)'

            % Update initial conditions for next period
            x0 = [OUT.S(end,:)'; OUT.D(end,:)'; OUT.N(end,:)'; OUT.Su(end,:)'];
            % Concatenate states and launch rates
            X = [X;x0'];
            LAM1 = [LAM1;lam(:,1)'];
            LAM2 = [LAM2;lam(:,2)'];
        end
        %% GMPHD block
        if strcmp(model_type, 'GMPHD')
            % Create solver guess
            solver_guess = lam(:,2)' + 100

            % Package all the GMPHD parameters into a single struct to pass to openAccessSolver
            GMPHD_all.t = t;
            GMPHD_all.GMPHD_params = GMPHD_params;
            GMPHD_all.shellstruct = shellstruct;
            GMPHD_all.GM_coll = GM_coll;
            GMPHD_all.i_collisions = i_collisions;
            GM1_last = GM1(end,:); %% Extract last element of GM1 to use in openAccessSolver

            % Solve for equilibrium launch rates
            lam(:,2) = openAccessSolver(solver_guess, launch_mask, 'GMPHD', GMPHD_all, GM1_last, "linear", econ_params, lam(:,1), n_workers)';

            % Update initial conditions for next period
            x0 = [Si';D_lt';D_st';D_lnt';Sui'];
            % Concatenate states and launch rates
            X = [X;x0'];
            LAM1 = [LAM1;lam(:,1)'];
            LAM2 = [LAM2;lam(:,2)'];
        end
    end

    elapsedTime = toc;
    fprintf('Time taken for period %d: %.2f seconds\n', timeIDX, elapsedTime);
end

%%%% WRITE OUT

if strcmp(model_type, 'GMPHD')
    %% Write outputs of model run as a CSV file
    S = X(:,1:VAR.N_shell); % Slotted objects
    Dlt = X(:,VAR.N_shell+1:2*VAR.N_shell); % Large trackable objects
    Dst = X(:,2*VAR.N_shell+1:3*VAR.N_shell); % Small trackable objects
    Dlnt = X(:,3*VAR.N_shell+1:4*VAR.N_shell); % LNTs
    Su = X(:,4*VAR.N_shell+1:5*VAR.N_shell); % Unslotted objects
    N_tot = S+Dlt+Dst+Dlnt+Su; % total population at each time instant for each shell

    %% Write out computed outputs to a folder with the same name as the stem of the input file
    dataTable = table(S', Dlt', Dst', Dlnt', Su', LAM1', LAM2', 'VariableNames', {'S', 'Dlt', 'Dst', 'Dlnt', 'Su', 'lSlotted', 'lUnslotted'});
    writetable(dataTable, ['scenarios/', stem, '/', filename]);
    fprintf('Wrote output to CSV');
end

%% MOCAT block
if strcmp(model_type, 'MOCAT')
    %% Write outputs of model run as a CSV file
    S = X(:,1:VAR.N_shell); % Slotted objects
    D = X(:,VAR.N_shell+1:2*VAR.N_shell); % Derelict objects
    N = X(:,2*VAR.N_shell+1:3*VAR.N_shell); % Debris objects
    Su = X(:,3*VAR.N_shell+1:4*VAR.N_shell); % Unslotted objects
    N_tot = S+D+N+Su; % total population at each time instant for each shell

    %% Write out computed outputs to a folder with the same name as the stem of the input file
    dataTable = table(S', D', N', Su', LAM1', LAM2', 'VariableNames', {'S', 'D', 'N', 'Su', 'lSlotted', 'lUnslotted'});
    writetable(dataTable, ['scenarios/', stem, '/', filename]);
    fprintf('Wrote output to CSV');

    %% Write out values necessary to compute collision rates in each shell. Loop over shells to fill a table in which rows are shells and columns are variables. The variables include all the inputs to the collision rate calculation in MOCAT4S.
    % Initialize the table
    varNames = {'Location', ...
        'SlottedRetire', 'SlottedDebris', 'SlottedDerelict', 'SlottedUnslotted', 'SlottedSlotted', ...
        'DerelictFromSlotted', 'DerelictDecay', 'DerelictSlotted', 'DerelictDebris', 'DerelictDerelict', 'DebrisSlotted', 'DerelictFromUnslotted', 'DerelictUnslotted', 'DebrisUnslotted', ...
        'UnslottedRetire', 'UnslottedDebris', 'UnslottedDerelict', 'UnslottedUnslotted', 'UnslottedSlotted'
        };
    numShells = VAR.N_shell;  % Replace N_shell with the actual number of shells if different
    parameterTable = array2table(zeros(0, numel(varNames)), 'VariableNames', varNames);

    % Define symbolic variables for populations in each shell k
    S = sym('S',[VAR.N_shell,1]);
    D = sym('D',[VAR.N_shell,1]);
    N = sym('N',[VAR.N_shell,1]);
    Su = sym('Su',[VAR.N_shell,1]);

    for k=1:VAR.N_shell
        % Flux within shell k
        rhok = densityexp(VAR.R02(k)); % no solar flux
        rvel_current_D=-rhok*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);
        rvel_current_N=-rhok*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);
        phi_SS = VAR.phi(1,1)/VAR.V(k);
        phi_SD = VAR.phi(1,2)/VAR.V(k);
        phi_SN = VAR.phi(1,3)/VAR.V(k);
        phi_SSu = VAR.phi(1,4)/VAR.V(k);
        phi_DD = VAR.phi(2,2)/VAR.V(k);
        phi_DN = VAR.phi(2,3)/VAR.V(k);
        phi_DSu = VAR.phi(2,4)/VAR.V(k);
        phi_NN = VAR.phi(3,3)/VAR.V(k);
        phi_NSu = VAR.phi(3,4)/VAR.V(k);
        phi_SuSu = VAR.phi(4,4)/VAR.V(k);
        % Parameters for slotted satellites in shell k
        A1 = -1/VAR.Dt; % A1 - decay rate
        A2 = -(VAR.delta+VAR.alpha)*phi_SN; % A2 - collisions with debris
        A3 = -(VAR.delta+VAR.alpha)*phi_SD; % A3 - collisions with derelicts
        A4 = -(1-VAR.slot_par)*VAR.alpha_active*phi_SS; % A4 - collisions between slotted satellites
        A5 = -VAR.alpha_active*phi_SSu; % A5 - collisions with unslotted satellites; needs to have factor of S/(Su+S) multiplied in later to recover time series of probability
        % Parameters for derelicts in shell k
        B1 = (1-VAR.P)/VAR.Dt;%S(k) % creation from retired slotted satellites
        B2 = rvel_current_D/VAR.Dhl;%D(k) % decay due to drag
        B3 = +VAR.delta*phi_SD;%D(k)*S(k) % collisions S-D
        B4 = -phi_DN;%(N(k))*D(k) % collisions D-N 
        B5 = -phi_DD;%(D(k))*D(k) % collisions D-D
        B6 = +VAR.delta*phi_SN;%N(k)*S(k) % collisions N-S
        B7 = (1-VAR.P)/VAR.Dt;%Su(k) % creation from retired unslotted satellites
        B8 = +VAR.delta*phi_DSu;%D(k)*Su(k) % collisions D-Su 
        B9 = +VAR.delta*phi_NSu;%N(k)*Su(k) % collisions N-Su
        % Parameters for unslotted satellites in shell k
        A1_Su = -1/VAR.Dt; % Decay rate of unslotted satellites
        A2_Su = -(VAR.delta+VAR.alpha)*phi_NSu; % Collisions with debris
        A3_Su = -(VAR.delta+VAR.alpha)*phi_DSu; % Collisions with derelicts
        A4_Su = -VAR.alpha_active*phi_SuSu; % Collisions between unslotted satellites
        A5_Su = -VAR.alpha_active*phi_SSu; % Collisions with slotted satellites; needs to have factor of Su/(Su+S) multiplied in later to recover time series of probability

        % Append the data for this shell to the table
        newRow = table(k, A1, A2, A3, A4, A5, ...
                    B1, B2, B3, B4, B5, B6, B7, B8, B9, ...
                    A1_Su, A2_Su, A3_Su, A4_Su, A5_Su, ...
                    'VariableNames', varNames);
        parameterTable = [parameterTable; newRow];
    end

    %% Make folder if it doesn't exist already
    mkdir(['scenarios/', stem]);

    %% Write parameterTable out to a CSV in the same folder as the output file
    writetable(parameterTable, ['scenarios/', stem, '/', stem, '--collision-parameters.csv']);

    % Write out VAR.K0 as a matrix into a CSV in the same folder as the output
    K0out = array2table(VAR.K0, 'VariableNames', {'S', 'Su', 'N', 'D'}, 'RowNames', {'S', 'Su', 'N', 'D'});
    K0out_filename = fullfile('scenarios', stem, 'fragmentation-parameters.csv');
    writetable(K0out, K0out_filename, 'WriteRowNames', true);
end
%% Write out parameter set with unique name
% Recreate a table to contain the parameters
parameter_table = [];

% Extract fields and values from VAR
fields_VAR = fieldnames(VAR);
values_VAR = struct2cell(VAR);
parameter_type_VAR = repmat({'MOCAT'}, numel(fields_VAR), 1);
parameter_table_VAR = table(parameter_type_VAR, fields_VAR, values_VAR, 'VariableNames', {'parameter_type', 'parameter_name', 'parameter_value'});

% Extract fields and values from econ_params
fields_econ = fieldnames(econ_params);
values_econ = struct2cell(econ_params);
parameter_type_econ = repmat({'econ'}, numel(fields_econ), 1);
parameter_table_econ = table(parameter_type_econ, fields_econ, values_econ, 'VariableNames', {'parameter_type', 'parameter_name', 'parameter_value'});

% Extract fields and values from GMPHD_params
fields_GMPHD = fieldnames(GMPHD_params);
values_GMPHD = struct2cell(GMPHD_params);
parameter_type_GMPHD = repmat({'GMPHD'}, numel(fields_GMPHD), 1);
parameter_table_GMPHD = table(parameter_type_GMPHD, fields_GMPHD, values_GMPHD, 'VariableNames', {'parameter_type', 'parameter_name', 'parameter_value'});

% Extract fields and values from constellation_params
fields_constellation = fieldnames(constellation_params);
values_constellation = struct2cell(constellation_params);
parameter_type_constellation = repmat({'constellation'}, numel(fields_constellation), 1);
parameter_table_constellation = table(parameter_type_constellation, fields_constellation, values_constellation, 'VariableNames', {'parameter_type', 'parameter_name', 'parameter_value'});

% Extract fields and values from launch_mask
fields_launch_mask = {'launch_mask'};
values_launch_mask = {launch_mask};
parameter_type_launch_mask = repmat({'launch_mask'}, numel(fields_launch_mask), 1);
parameter_table_launch_mask = table(parameter_type_launch_mask, fields_launch_mask, values_launch_mask, 'VariableNames', {'parameter_type', 'parameter_name', 'parameter_value'});

% Concatenate all the tables
scenario_table = [parameter_table_VAR; parameter_table_GMPHD; parameter_table_econ; parameter_table_constellation; parameter_table_launch_mask];

% Write the unique name to a text CSV file with the unique name as filename and scenario table as contents. The file is saved in the "scenarios" folder under a sub-folder with the same unique name. The sub-folder must be created first.
writetable(scenario_table, ['scenarios/', stem, '/', stem, '--parameters.csv']);

big_timer_elapsed = toc(big_timer); % Stop the timer and get the elapsed time
fprintf('Total time taken: %.2f seconds\n', big_timer_elapsed); % Print the time taken
