% This is a solver for the following problem:
% * Given one starting point q0 and one end point qn
% * Deliver a package from q0 to qn
% * There are two drones nearby at position p0 and p1, one with speed v0 
% and one with speed v1 with v0 < v1
% * Assume that drone 0 can get to q0 faster than drone 1
% * After some time, drone 0 will hand the package to drone 1 at some point
% q1 (q1 is the only variable in this problem)
% * t0: time for drone 0 to fly from p0 to q0, a constant
% * t01: time for drone 0 to fly from q0 to q1 = time for dron 1 to fly
% from p1 to q1
% * t12: time for drone 1 to fly from q1 to qn
% * Minimize t01 + t12
% * st. ||q1 - q0||/v0 - ||q1 - p1||/v1 + t0 = 0 
% * t01 = ||q1 - q0||/v0 = ||q1 - p1||/v1 - t0
% * t12 = ||qn - q1||/v1
syms xq1 yq1; %q1
syms xq0 yq0; %q0
syms xp1 yp1; %p1
syms xqn yqn; %qn
syms v0 v1; %velocity
syms t0; %constant
syms lamb; %langrange mult
denom_q1p1 = sqrt((xq1 - xp1)^2 + (yq1 - yp1)^2);
denom_q1qn = sqrt((xq1 - xqn)^2 + (yq1 - yqn)^2);
denom_q1q0 = sqrt((xq1 - xq0)^2 + (yq1 - yq0)^2);

eqnx = (1-lamb)*(xq1 - xp1)/(v1*denom_q1p1) + (xq1-xqn)/(v1*denom_q1qn) +...
    lamb*(xq1-xq0)/(v0*denom_q1q0)==0;
eqny = (1-lamb)*(yq1 - yp1)/(v1*denom_q1p1) + (yq1-yqn)/(v1*denom_q1qn) +...
    lamb*(yq1-yq0)/(v0*denom_q1q0)==0;
constraint = denom_q1q0/v0 - denom_q1p1/v1 + t0 == 0;

[solx, soly] = solve([eqnx, eqny, constraint], [xq1, yq1, lamb]);