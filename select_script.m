function output = select_script(run)

% This file contains Matlab codes for the stability and sensitivity
% analysis of channel flows, Poiseuille and Couette flows. These
% are parts of the tutorial: "Analysis of fluid systems: stability,
% receptivity, sensitivity" by Peter Schmid and Luca Brandt,
% published in Applied Mechanics Reviews, 66(2), 2014.
%
% The main programs are
%
% TransientGrowth.m    : compute the transient growth curve G(t)
% OptimalDisturbance.m : compute the the optimal initial condition 
%                          and the corresponding flow response
% Neutral_a_Re.m       : compute the maximum optimal growth and the
%                          least stable eigenvalue in the Reynolds-
%                          alpha plane
% Neutral_alpha_beta.m : compute the maximum optimal growth and the
%                          least stable eigenvalue in the alpha-
%                          beta plane
% Resolvent.m          : compute the resolvent norm for real and
%                          complex frequency omega
% NumRange.m           : compute the spectrum and numerical range
%                          of the stability operator
%
% To execute, replace the argument 'run' by the string obtained
% from the corresponding function name, sans the suffix.

switch run
  case 'TransientGrowth'
    output = TransientGrowth();
  case 'OptimalDisturbance'
    output = OptimalDisturbance();
  case 'Neutral_a_Re'
    output = Neutral_a_Re();
  case 'Neutral_alpha_beta'
    output = Neutral_alpha_beta();
  case 'Resolvent'
    output = Resolvent();
  case 'NumRange'
    output = NumRange();
  otherwise
    warning('AMR:nodemo', 'No such demo');
end

end