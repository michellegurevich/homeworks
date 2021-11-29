%%% Research Skills, Matlab final
%%% Michelle Gurevich

% final_gurevich.m
% Submission for final project for Matlab portion of course

% Author: Michelle Gurevich
% Date: 4 June 2021

clear
clc
close all

%% Part 1: The lifetime of atoms in a magneto-optical-trap (MOT)

% Import data file and plot contents
Data = importdata( 'MOT_Lifetime.txt' );
figure
plot(Data(:,1), Data(:,2))
grid on
xlabel( 'Time / seconds' );
ylabel( 'Number of atoms' );
title( 'Lifetime of atoms in a magneto-optical-trap' );

% Plot number of atoms over time after decay begins
figure
index = find(Data(:,1) == 0,2);
NDecay = Data(index:end,:);
plot(NDecay)
grid on
xlabel( 'Time / seconds' );
ylabel( 'Number of atoms' );
title( 'Lifetime of atoms in MOT after decay begins' );

% Use ExpDecay function to fit data
myfit = fittype( 'ExpDecay(t, tau, A, B)', 'independent', 't', ...
    'coefficients',  { 'tau', 'A', 'B' });
fobj = fit( NDecay(:,1), NDecay(:,2), myfit, 'StartPoint', [2.5 10e5 0]);
figure
% Plot data with ExpDecay fit
plot( fobj, NDecay(:,1), NDecay(:,2) )
grid on
xlabel( 'Time / seconds' );
ylabel( 'Number of atoms' );
title( 'ExpDecay fit of atoms in MOT with decay' );

%% Part 2: Dynamics of a Particle in a Double-Well Potential

% Task One
% construct ladder operators
D = 25; %truncate dimension for numerical impl
a = diag(sqrt(1:D-1), 1);
adagger = a';
    
% construct position and momentum operators
hbar = 1;
m = 0.1;
omega = 1;
p = -1i * sqrt((hbar*m*omega)/2) * (a-adagger);
x = sqrt(hbar/(2*m*omega)) * (a+adagger);
    
% Task Two
% build hamiltonian from eq 2
H = p^2 / (2*m) + (1/2)* m * omega^2 * x^2;
% build hamiltonian from eq 5
Hhat = hbar * omega * (adagger * a + 1/2 * eye(D));
    
% the position and momentum operators use the states preceding and
% following and so by truncating at a finite value we miss out on a
% contribution to the value of the value at truncation index

% Task Three
[basis, D] = eig(Hhat);

% the columns of basis correspond to |psi> the eigenstates of H hat and the
% eigenenergies n are the diagonal values of D, where each eigenenergie is
% an eigenvalue of the system and is given by hbar * omega (n + 1/2)

x1 = linspace(-10,10,200);
y1 = ones(1,200);
V = 1/2 * m * omega^2 * x1.^2 % compute the potential of the harmonic osc

% plot the harmonic potential and and the five lowest energy levels
figure
plot(x1, V, x1, D(1,1)*y1, x1, D(2,2)*y1, x1, D(3,3)*y1, x1, D(4,4)*y1, ... 
    x1, D(5,5)*y1)
xlabel( 'x' );
ylabel( 'Energy levels' );
title( 'Harmonic potential and lowest five energy levels' );
legend( 'Harmonic potential', 'n=1', 'n=2', 'n=3', 'n=4', 'n=5')