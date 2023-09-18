function [econParams] = set_econ_parameters(VAR)
%%%%%
% Function to generate economics parameters of integrated assessment model
%%%%%

%% Lifetime of a satellite [years]. Stored here to simplify excess returns calls; should be checked to ensure consistency with physical model          
% Default values:
% - MOCAT4S: 5 years
econParams.satLifetime = 5;

%% Regulated disposal rule time [years], default is 5
econParams.disposalTime = 5;

%% Discount rate [pp/year]
econParams.discountRate = 0.05;

%% Parameters specific to linear revenue model [$/year] . A satellite facing no competition earns 750,000 $/year in revenues.
econParams.intercept = 7.5e5; % [$/year]
econParams.coef = 1.0e2; % [$/satellite/year]
econParams.tax = 0.0; % [%] tax rate on shell-specific pc

%% Cost for a single satellite to use any shell [$] .

%%% Cost of delta-v [$/km/s]
econParams.delta_v_cost = 1000;

%%% Price of lift [$/kg]  
% Default is $5000/kg based on price index calculations in Corrado, Cropper, Rao (2023)
econParams.lift_price = 5000;
