%% Stochastic simulation of the dynamics of neurofilaments through 
% axonal internodes and nodes of Ranvier
% Calls functions: initialize.m, content.m, contenton,m, update.m

clear all
close all

va = 0.5; vr = -0.5;

% Kinetic rates from previous studies
g01  = 0.064;
g10  = 0.14;
goff = 0.0045;
gon0 = 0.0001;   % is actually modulated
gar  = 0.0000042;
gra  = 0.000014;

% Length axon
axmin = 0;
axmax = 800; 

% Here could prescribe a node location, or go beyond axmax for an internode
% simulation
n1 = 805;
n2 = 805;

% Time parameters
dt      = 1;           % time step (seconds)
t_end   = 25*60;       % end time (seconds)
t_begin = 0;           % start time
delta_t = 60;          % time step (seconds)

number_mt    = 100;    % number of microtubules
seats_per_mt = 5;      % lanes of traffic per microtubule
nrunners     = 213500; % number of neurofilaments


% NF length distribution
nf_leng   = 1;           % 0 if uniform, 1 if exponential, 2 if discrete data distribution
maxlength = 6;           % max length of NFs: for uniform length
avgleng   = maxlength;   % average of exponential distribution for NF lengths
if nf_leng == 0
    leng  = maxlength.*ones(1,nrunners);
elseif nf_leng == 1
    leng  = ceil(exprnd(avgleng,1,nrunners));
else                     % given an empirical NF length distribution
    table = csvread('FilamentLengthDistribution.csv');
    lengs = table(:,1);  % binned 1-micron NF lengths
    freqs = table(:,2);  
    lengs = lengs + 1;   % no 0 lengths

    freqs(end) = 1-sum(freqs(1:end-1));
    x0    = randsmpl(freqs, 1, nrunners); % draw from discrete distribution
    leng  = lengs(x0);
    
end

nr_times  = floor(t_end/delta_t);
xtable    = zeros(nr_times,nrunners);
volume_table       = zeros(nr_times,axmax);
occupancy_table    = zeros(nr_times,axmax);
nstate_table       = zeros(nr_times,nrunners);
nrunners_table     = zeros(1,nr_times);
gon_vect_table     = zeros(nr_times,nrunners);
var_vect_table     = zeros(nr_times,nrunners);

% Initialize NF center locations and kinetic states, and calculate NF
% content in space
[x,nstate]  = initialize(g10,g01,goff,gon0,gar,gra,nrunners,axmin,axmax); 
[volume]    = content(x,leng,axmax,nrunners); 
[occupancy] = contenton(x,leng,axmax,nrunners,nstate); 

%% Time loop

t = t_begin;
n = 0;
volume_taver = zeros(1,axmax);

while t<t_end
    
    t/delta_t
    
    [xnew,nstatenew,tnew,volumenew,occupancynew,gon_vect,var_vect] = update_nfs(x,leng,nstate,...
        volume,occupancy,t,t+delta_t,nrunners,number_mt,seats_per_mt,...
        goff,g10,g01,gar,gra,va,vr,axmax,axmin,n1,n2,dt);
    
    t = tnew;
    x = xnew;
    nstate    = nstatenew;
    volume    = volumenew;
    occupancy = occupancynew;
    
    x_current = x;
    n = n+1; 
    
    nstate_table(n,:)    = nstatenew;
    gon_vect_table(n,:)  = gon_vect;
    var_vect_table(n,:)  = var_vect;
    xtable(n,:)          = x_current;
    
    volume_table(n,:)    = volume;
    occupancy_table(n,:) = occupancy;
    nrunners_table(1,n)  = nrunners;
    
end

save('node_test.mat')

%% Plotting
for i = 1:size(xtable,1)
    plot(volume_table(i,:),'Linewidth',2)
    hold on
    plot(occupancy_table(i,:),'Linewidth',2,'Color',[0 0.5 0]) 
%     xlim([0 axmax+5]);
    xlim([350 450]);
    ylim([0 3000]);
    
    % Plot vertical lines if there is a node
    plot([n1 n1], ylim,'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5]) % Plot Vertical Line
    plot([n2 n2], ylim,'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5]) % Plot Vertical Line
    hold off
    
    % Labeling
    xlabel('Distance along axon (\mu m)')
    ylabel('NF content')
    title(['Minute ' num2str(i)]) %
    legend('All NFs','On-track NFs')%'All NFs','Time-averaged NFs'
    set(gca,'FontSize',16)
    drawnow
    
    pause(0.1)
end
