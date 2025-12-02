% Li-ion lumped model
% Runs open-circuit and simple kinetic+ohmic model vs SoC

% change xnomin/max and xprmin/max to match pristine cell data
% adjust LLI value to match degraded cell data 

clearvars
close all
clc

%% ----------------- USER / PHYSICAL CONSTANTS -----------------
T    = 303;              % K
R    = 8.314;            % J/(mol K)
F    = 96485;            % C/mol
alfa = 0.5;              % charge transfer symmetry factor (dimensionless)
b    = R*T/F/alfa;       % prefactor in overpotential expression (V)

%% ----------------- DEGRADATION / CAPACITY -----------------
LLI    = 0.268;           % Loss of lithium inventory (fraction), 0..1
C_n    = 41 / 1000 * (1 - LLI);% Nominal battery capacity [Ah] (affected by LLI)

%% ----------------- NEGATIVE ELECTRODE PARAMETERS --------------
xnomin = 1e-3;           % stoichiometry lower bound (neg)
xnomax = 0.8;            % stoichiometry upper bound (neg)
nn     = 2;              % exponent to map stoichiometry->activity
Fkn    = 10;             % scale for j0 (A) - empirical parameter
E0n    = 0.1;            % reference potential negative electrode [V]

%% ----------------- POSITIVE ELECTRODE PARAMETERS --------------
xprmin = 0.03;           % positive electrode reference stoich min
xprmax = 0.9;
xpomin = 1 - xprmax;     % mapping used by original code
xpomax = 1 - xprmin;
np     = 5;              % exponent to map stoichiometry->activity
Fkp    = 100;            % scale for j0 (A) - empirical parameter
E0p    = 3.95;           % reference potential positive electrode [V]

%% ----------------- ELECTROLYTE / OHMICS ----------------------
Re     = 15e-3;          % electrolyte/pack resistance [Ohm]
% You may increase Re with LLI to mimic SEI growth increasing internal resistance:
% Re = Re*(1 + k * LLI) where k is chosen empirically.

%% ----------------- SIMULATION GRID / CURRENT ------------------
SoC    = 0:0.001:1;       % state-of-charge vector (0..1)
npts   = numel(SoC);

% Discharge/charge choice
DC     = 1;              % = 1 discharge, = -1 charge
Crate  = 2;              % C-rate [1/h] (2C means discharge in 0.5 h)

% Current magnitude (A).
% NOTE: C_n is [Ah], Crate is [1/h] -> C_n*Crate is [A].
j      = DC * C_n * Crate;   % signed current (A). positive => discharge here.

%% ----------------- ACTIVITIES (simple power-law mapping) ------------
% Activities computed from stoichiometry; powers nn/np used to mimic
% concentration dependence (simple phenomenological model).
ano = ((xnomax - xnomin) * (1 - LLI) .* SoC + xnomin) .^ nn;       % Li in negative (active)
anr = (1 - xnomin - (xnomax - xnomin) * (1 - LLI) .* SoC) .^ nn;   % vacant in neg
apo = ((xpomax - xpomin) * (1 - SoC) * (1 - LLI) + xpomin) .^ np;  % Li in positive
apr = (1 - xpomin - (xpomax - xpomin) * (1 - SoC) * (1 - LLI)) .^ np;% vacant in pos

% Numerical safety: avoid log(0) or extremely small values
tiny = 1e-12;
ano = max(ano, tiny);
anr = max(anr, tiny);
apo = max(apo, tiny);
apr = max(apr, tiny);

%% ----------------- OCV calculation (simple Nernst-like) ------------
Ep = E0p + (R*T/F) * log(apr ./ apo);   % positive electrode potential vs Li/Li+
En = E0n + (R*T/F) * log(anr ./ ano);   % negative electrode potential vs Li/Li+
OCV = Ep - En;                          % open circuit cell voltage

%% ----------------- Kinetics: exchange-like currents j0 (A) ----------
% j0 ~ Fk * a_ox^(1-alfa) * a_red^alfa  (very simplified)
jp = Fkp * apo.^(1 - alfa) .* apr.^alfa;  % pos electrode effective j0 (A)
jn = Fkn * ano.^(1 - alfa) .* anr.^alfa;  % neg electrode effective j0 (A)

% Avoid jp/jn being zero
jp = max(jp, tiny);
jn = max(jn, tiny);

%% ----------------- OVERPOTENTIALS --------------------------------
% Use signed current in asinh so that charge/discharge flip sign correctly.
% Butler-Volmer small/large current approx: eta = ± b * asinh(j/(2*j0))
etap = - b * asinh( j ./ (2 * jp) );   % positive electrode overpotential
etan =   b * asinh( j ./ (2 * jn) );   % negative electrode overpotential

% Ohmic drop (electrolyte + internal resistance), sign follows current.
etae = Re * j * ones(size(SoC));

% Terminal voltage (accounting for electrode overpotentials and ohmic drop)
V_terminal = (Ep + etap) - (En + etan) - etae;

%% ----------------- DIAGNOSTICS: check for NaNs/Inf ----------------
if any(~isfinite([Ep, En, OCV, etap, etan, etae, V_terminal]))
    warning('Some computed values are non-finite. Check stoichiometry / activity ranges and tiny clamping.');
end

%% ----------------- PLOTS ------------------------------------------
mAh_exchanged = (1 - SoC) * C_n * 1000; % mAh on x-axis (like original)

figure('Name', 'Open-circuit condition'); hold on
plot(mAh_exchanged, Ep,  'LineWidth', 2);
plot(mAh_exchanged, En,  'LineWidth', 2);
plot(mAh_exchanged, OCV, 'LineWidth', 2);
legend('Pos OCP','Neg OCP','OCV','Location','Best');
grid on
xlabel('Exchanged charge [mAh]');
ylabel('Voltage [V]');
title(sprintf('OCV curves (LLI=%.3f, Crate=%.2fC)', LLI, Crate))

figure('Name', 'Overpotentials'); hold on
plot(mAh_exchanged, etap, 'LineWidth', 2);
plot(mAh_exchanged, etan, 'LineWidth', 2);
plot(mAh_exchanged, etae, 'LineWidth', 2);
legend('Pos \eta','Neg \eta','Electrolyte \eta','Location','Best');
grid on
xlabel('Exchanged charge [mAh]');
ylabel('Voltage [V]');
title('Overpotentials')

figure('Name', 'Battery operation'); hold on
plot(mAh_exchanged, Ep + etap, 'LineWidth', 2);
plot(mAh_exchanged, En + etan, 'LineWidth', 2);
plot(mAh_exchanged, V_terminal, 'LineWidth', 2);
legend('Pos Potential','Neg Potential','Terminal V','Location','Best');
grid on
xlabel('Exchanged charge [mAh]');
ylabel('Voltage [V]');
title('Potentials and Terminal Voltage')

figure('Name', 'Activities'); hold on
plot(SoC, apo, 'LineWidth', 2);
plot(SoC, apr, 'LineWidth', 2);
plot(SoC, ano, 'LineWidth', 2);
plot(SoC, anr, 'LineWidth', 2);
legend('apo','apr','ano','anr','Location','Best');
grid on
xlabel('SoC [-]');
ylabel('Activities [-]');
title('Activity profiles')

%% ----------------- PRINT SUMMARY ---------------------------------
fprintf('Summary:\n');
fprintf(' C_n = %.3f Ah, Crate = %.2fC -> |I| = %.3f A (signed I = %.3f A)\n', C_n, Crate, abs(C_n*Crate), j);
fprintf(' Min terminal V = %.3f V, Max terminal V = %.3f V\n', min(V_terminal), max(V_terminal));
