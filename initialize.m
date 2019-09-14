%% Initialize NF center localions and kinetic states

function [x,nstate] = initialize(g10,g01,goff,gon0,gar,gra,nrunners,...
                                 axmin,axmax)

q1  = g10/g01;
q2  = goff/gon0;
q3  = gar/gra;
pa  = 1/((1+q1*(1+q2))*(1+gar/gra)); % state 1: anterograde moving
pa0 = pa*q1; % state 2: anterograde pausing on track
pr  = pa*q3; % state 3: anterograde pausing off track
pr0 = pr*q1; % state 4: retrograde moving
paoff = pa*q1*q2; % state 5: retrograde pausing on track
% proff = pr*q1*q2*q3; % state 6: retrograde pausing off track

nstate = zeros(1,nrunners);
x      = zeros(1,nrunners);

r1 = rand(1,nrunners);
r2 = rand(1,nrunners);

% Generate uniformly-spaced NF locations at start
% Assign states according to stationary distribution
for i = 1:nrunners

    x(i) = axmin + r1(i)*(axmax-axmin); 
    
    if r2(i)<=pa
       nstate(i) = 1;
    elseif r2(i)<=pa+pa0
        nstate(i) = 2;
    elseif r2(i)<=pa+pa0+paoff
        nstate(i) = 3;
    elseif r2(i)<=pa+pa0+paoff+pr
        nstate(i) = 4;
    elseif r2(i)<=pa+pa0+paoff+pr+pr0
        nstate(i) = 5;
    else
        nstate(i) = 6;
    end
    
end



