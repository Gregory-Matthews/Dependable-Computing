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
p_o_sol(t) = S.p_2 + S.p_1;

% p_o_sol_numeric = subs(p_o_sol, [lambda,mu, c], [0.001, 1/24, 1]);
% ezplot(p_o_sol_numeric, [0 time]);
% hold on;
% p_o_sol_numeric = subs(p_o_sol, [lambda,mu, c], [0.001, 1/48, 1]);
% ezplot(p_o_sol_numeric, [0 time]);
% hold on;
% p_o_sol_numeric = subs(p_o_sol, [lambda,mu, c], [0.001, 1/96, 1]);
% ezplot(p_o_sol_numeric, [0 time]);

% grid on;
% legend('\mu = 1/24', '\mu = 1/48', '\mu = 1/96')
% xlabel('Time (t)'); ylabel('Availability, A(t)');
% title('Dual-Modular Redundant System Availability');

%% Solve DMR model Reliability
% p_2: probability that both modules are operating
% p_1: probability that one of the two modules are operating
% p_f: probability that the system has failed

syms p_2(t) p_1(t) p_f(t)
ode1 = diff(p_2) == 2*mu*p_1 - 2*lambda*p_2;
ode2 = diff(p_1) == (c+1)*lambda*p_2 - (2*lambda + 2*mu)*p_1;
ode3 = diff(p_f) == (1-c)*lambda*p_2 + 2*lambda*p_1;
odes = [ode1; ode2; ode3];

cond1 = p_2(0) == 1; cond2 = p_1(0) == 0; cond3 = p_f(0) == 0;
conds = [cond1; cond2; cond3];

S = dsolve(odes, conds);
% p_2_sol(t) = S.p_2; p_1_sol(t) = S.p_1; p_f_sol(t) = S.p_f;
% p_o_sol(t) = S.p_2 + S.p_1;
% p_o_sol_numeric = subs(p_o_sol, [lambda,mu, c], [0.001, 1/48, 0.8]);
% ezplot(p_o_sol_numeric, [0 time]);
% hold on;
% p_o_sol_numeric = subs(p_o_sol, [lambda,mu, c], [0.001, 1/48, 0.9]);
% ezplot(p_o_sol_numeric, [0 time]);
% hold on;
% p_o_sol_numeric = subs(p_o_sol, [lambda,mu, c], [0.001, 1/48, 0.95]);
% ezplot(p_o_sol_numeric, [0 time]);

% grid on;
% legend('C = 0.8', 'C = 0.9', 'C = 0.95')
% xlabel('Time (t)'); ylabel('Reliability, R(t)');
% title('Dual-Modular Redundant System Reliability');


%% Solve TDTMR vs. QUAD Reliability

% TDTMR Analysis
syms p_3(t) p_2(t) p_1(t) p_ud(t) p_f(t)
ode1 = diff(p_3) == -3*lambda*p_3;
ode2 = diff(p_2) == 3*lambda*c*p_3 -2*lambda*p_2;
ode3 = diff(p_1) == 2*lambda*p_2 - lambda*p_1;
ode4 = diff(p_ud) == 3*lambda*(1-c)*p_3 - 2*lambda*p_ud;
ode5 = diff(p_f) == 2*lambda*(1-c)*p_2 + lambda*p_1 + 2*lambda*p_ud;
odes = [ode1; ode2; ode3; ode4; ode5];

cond1 = p_3(0) == 1; cond2 = p_2(0) == 0; cond3 = p_1(0) == 0;
cond4 = p_ud(0) == 0; cond5 = p_f(0) == 0;
conds = [cond1; cond2; cond3; cond4; cond5];

S = dsolve(odes, conds);
p_3_sol(t) = S.p_3;p_2_sol(t) = S.p_2;p_1_sol(t) = S.p_1; 
p_ud_sol(t) = S.p_ud; p_f_sol(t) = S.p_f;
R_sol(t) = S.p_3 + S.p_2 + S.p_1 + S.p_ud;

% R_sol_numeric = subs(R_sol, [lambda, c], [0.0004, 0.98]);
% ezplot(R_sol_numeric, [0 time]);
% hold on;
% R_sol_numeric = subs(R_sol, [lambda, c], [0.005, 0.9]);
% ezplot(R_sol_numeric, [0 time]);
% hold on;
% R_sol_numeric = subs(R_sol, [lambda, c], [0.0025, 0.9]);
% ezplot(R_sol_numeric, [0 time]);
% hold on;

% QUAD Analysis
syms p_4(t) p_3(t) p_2(t) p_1(t) p_ud3(t) p_ud2(t) p_f(t)
ode1 = diff(p_4) == -4*lambda*p_4;
ode2 = diff(p_3) == 4*lambda*c*p_4 - 3*lambda*p_3;
ode3 = diff(p_2) == 3*lambda*c*p_3 - 2*lambda*p_2;
ode4 = diff(p_1) == 2*lambda*c*p_2 - lambda*p_1;
ode5 = diff(p_ud3) == 4*lambda*(1-c)*p_4 - 3*lambda*p_ud3;
ode6 = diff(p_ud2) == 3*lambda*(1-c)*p_3 + 3*lambda*c*p_ud3 - 2*lambda*p_ud2;
ode7 = diff(p_f) == 3*lambda*(1-c)*p_ud3 + 2*lambda*(1-c)*p_2 + lambda*p_1 + 2*lambda*p_ud2;
odes = [ode1; ode2; ode3; ode4; ode5; ode6; ode7];

cond1 = p_4(0) == 1; cond2 = p_3(0) == 0; cond3 = p_2(0) == 0; 
cond4 = p_1(0) == 0; cond5 = p_ud3(0) == 0; 
cond6 = p_ud2(0) == 0; cond7 = p_f(0) == 0;
conds = [cond1; cond2; cond3; cond4; cond5; cond6; cond7];

S = dsolve(odes, conds);
p_4_sol(t) = S.p_4; p_3_sol(t) = S.p_3;p_2_sol(t) = S.p_2;
p_1_sol(t) = S.p_1; p_ud3_sol(t) = S.p_ud3; p_ud2_sol(t) = S.p_ud2; 
p_f_sol(t) = S.p_f;
R_sol(t) = S.p_4 + S.p_3 + S.p_2 + S.p_1 + S.p_ud3 + S.p_ud2;

% R_sol_numeric = subs(R_sol, [lambda, c], [0.001, 0.96]);
% ezplot(R_sol_numeric, [0 time]);
% hold on;
% R_sol_numeric = subs(R_sol, [lambda, c], [0.005, 0.9]);
% ezplot(R_sol_numeric, [0 time]);
% hold on;
% R_sol_numeric = subs(R_sol, [lambda, c], [0.0025, 0.9]);
% ezplot(R_sol_numeric, [0 time]);

% grid on;
% xlabel('Time (t)'); ylabel('Reliability, R(t)');
% 
% yt=get(gca,'YTick');
% ylab=num2str(yt(:), '%15.8f');
% set(gca,'YTickLabel',ylab); 
% 
% legend('TDTMR: \lambda = 0.001', 'TDTMR: \lambda = 0.005', 'TDTMR:\lambda = 0.0025','QUAD: \lambda = 0.001', 'QUAD: \lambda = 0.005', 'QUAD: \lambda = 0.0025')
% title('TDTMR System Reliability');
