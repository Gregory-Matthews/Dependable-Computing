% ECEC 520 - Assignment 1
% Author: Naga Kandasamy, Date: 4/12/17
% Edits By: Greg Matthews, Date: 4/24/17
%%
clear all; close all; clc;
syms lambda mu c; % The failure rate associated with a module
time = 4; % Number of hours of system operation
 
%% Solve DMR model Availability
% p_2: probability that both modules are operating
% p_1: probability that one of the two modules are operating
% p_f: probability that the system has failed

syms p_2(t) p_1(t) p_f(t)
ode1 = diff(p_2) == 2*mu*p_1 - 2*lambda*p_2;
ode2 = diff(p_1) == (c+1)*lambda*p_2 - (2*lambda + 2*mu)*p_1 + 2*mu*p_f;
ode3 = diff(p_f) == (1-c)*lambda*p_2 + 2*lambda*p_1 - 2*mu*p_f;
odes = [ode1; ode2; ode3];

cond1 = p_2(0) == 1; cond2 = p_1(0) == 0; cond3 = p_f(0) == 0;
conds = [cond1; cond2; cond3];

S = dsolve(odes, conds);
p_2_sol(t) = S.p_2; p_1_sol(t) = S.p_1; p_f_sol(t) = S.p_f;
p_o_sol(t) = S.p_2 + S.p_1
