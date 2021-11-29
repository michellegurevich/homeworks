%%% Research Skills, ExpDecay file
%%% Michelle Gurevich

% ExpDecay.m
% Supporting file for Part 1 of Matlab final

% Author: Michelle Gurevich
% Date: 4 June 2021

%% ExpDecay function

function fit_decay = ExpDecay(t, tau, A, B)
    fit_decay = A * exp(-t/tau) + B;
end
