function [xnew,nstatenew,tnew,volumenew,occupancynew,gon_vect,var_vect] = update_nfs(x,leng,nstate,volume,occupancy,t1,t2,nrunners,...
    number_mt,seats_per_mt,goff,g10,g01,gar,gra,va,vr,axmax,axmin,n1,n2,dt)

temp  = t1;
xtemp = x;
nstatetemp = nstate;

% Diffusion/effective cross-sectional area parameters
D = 8*10^(-6); 
anf_inter  = 4.2*10^(-3);  amt_inter = 2.91*10^(-2); % internode
anf_node   = 1.74*10^(-3); amt_node  = 1.19*10^(-2); % node


while temp<=t2
    
    r1 = rand(1,nrunners);
    
    for i = 1:nrunners
        
        n  = round(x(i));                 % bin of nf_i
        
        if (x(i)+leng(i)/2 - (floor(x(i)+leng(i)/2))) >= 0.5
            nl = round(x(i)-leng(i)/2);   % bin of the left end of nf_i
            nr = round(x(i)+leng(i)/2)+1; % bin of the right end of nf_i
        else
            nl = round(x(i)-leng(i)/2)-1; % bin of the left end of nf_i
            nr = round(x(i)+leng(i)/2);   % bin of the right end of nf_i
        end

        if nr>axmax
            nr = nr-axmax;
        elseif nr<=axmin
            nr = nr+axmax;
        end
        
        if n>axmax
            n = n-axmax;
        elseif n<=axmin
            n = n+axmax;
        end
        
        if nl<=axmin
            nl = nl+axmax;
        elseif nl>axmax
            nl = nl-axmax;
        end

            modul_loc = x(i);  % gamma_on modulation at the center of the NF
        
        if number_mt*seats_per_mt-occupancy(n)>=0
            if modul_loc>=n1 && modul_loc<=n2  %x(i)>=n1 && x(i)<=n2
                    gon = (4*D/seats_per_mt)*((number_mt*seats_per_mt-occupancy(n))/(number_mt*amt_node+volume(n)*anf_node));

            else

                    gon = (4*D/seats_per_mt)*((number_mt*seats_per_mt-occupancy(n))/(number_mt*amt_inter+volume(n)*anf_inter));

            end
        else  % if no available seats
            gon = 0;
        end
        
        gon_vect(i) = gon;
        
        % nr
        if number_mt*seats_per_mt-occupancy(nr)>=0
            va0 = va*(1-occupancy(nr)/(number_mt*seats_per_mt));
        else  % if no available seats
            va0 = 0;
        end
        
        % nl
        if number_mt*seats_per_mt-occupancy(nl)>=0
            vr0 = vr*(1-occupancy(nl)/(number_mt*seats_per_mt));
        else   % if no available seats
            vr0 = 0;
        end
        
        if nstate(i) == 1
            var_vect(i) = va0;
        elseif nstate(i) == 4
            var_vect(i) = vr0;
        else
            var_vect(i) = 0;
        end
        
        volume0    = volume;
        occupancy0 = occupancy;
        
        if nstate(i) == 1
            if r1(i)<=g10*dt
                nstatetemp(i) = 2;
                % 1->2, no change in volume or occupancy
            else
                xtemp(i) = x(i)+va0*dt;               % anterograde movement
                % 1->1
                if (x(i)+leng(i)/2 - (floor(x(i)+leng(i)/2))) >= 0.5
                    k1 = (round(xtemp(i)-leng(i)/2));
                else
                    k1 = (round(xtemp(i)-leng(i)/2)-1);
                end
                
                if nl~=k1   % moved by at least a bin, to the right                 
                    if nl+1>axmax
                        nlp = nl-axmax;
                    else
                        nlp = nl;
                    end
                    
                    volume0(nlp+1)    = volume0(nlp+1) - 1;
                    volume0(nr)       = volume0(nr) + 1;
                    occupancy0(nlp+1) = occupancy0(nlp+1) - 1;
                    occupancy0(nr)    = occupancy0(nr) + 1;
                end
            end
        elseif nstate(i)==2
            if r1(i)<=g01*dt
                nstatetemp(i)=1;
                x(i)=x(i)+va0*dt;
                % 2 -> 1, switches to moving anterograde
                if (x(i)+leng(i)/2 - (floor(x(i)+leng(i)/2))) >= 0.5
                    k1 = (round(xtemp(i)-leng(i)/2));
                else
                    k1 = (round(xtemp(i)-leng(i)/2)-1);
                end
                if nl~=k1   % moved by at least a bin, to the right                    
                    if nl+1>axmax
                        nlp = nl-axmax;
                    else
                        nlp = nl;
                    end
                    
                    volume0(nlp+1)    = volume0(nlp+1) - 1;
                    volume0(nr)       = volume0(nr) + 1;
                    occupancy0(nlp+1) = occupancy0(nlp+1) - 1;
                    occupancy0(nr)    = occupancy0(nr) + 1;
                end
            elseif r1(i)<=(g01+goff)*dt
                nstatetemp(i) = 3;
                % 2-> 3, one less in the on-track category                
                if nl+1>axmax
                    nlp = nl-axmax;
                else
                    nlp = nl;
                end
                if nr-1<=axmin
                    nrp = nr+axmax;
                else
                    nrp = nr;
                end
                
                for k=nlp+1:nrp-1 % Subtract from the old locations
                    if k>axmax
                        kk = k-axmax;
                    elseif k<1
                        kk = k+axmax;
                    else
                        kk=k;
                    end
                    occupancy0(kk) = occupancy0(kk)-1;
                end
                
            elseif r1(i)<=(g01+goff+gar)*dt
                nstatetemp(i) = 5;
                % 2 -> 5, no change
            end
        elseif nstate(i)==3
            if r1(i)<=gon*dt
                nstatetemp(i) = 2;
                % 3 -> 2, one more in the on-track category
                if nl+1>axmax
                    nlp = nl-axmax;
                else
                    nlp = nl;
                end
                if nr-1<=axmin
                    nrp = nr+axmax;
                else
                    nrp = nr;
                end
                
                for k=nlp+1:nrp-1 % Add to the locations
                    if k>axmax
                        kk = k-axmax;
                    elseif k<1
                        kk = k+axmax;
                    else
                        kk=k;
                    end
                    occupancy0(kk) = occupancy0(kk)+1;
                end
            elseif r1(i)<=(gon+gar)*dt
                nstatetemp(i) = 6;
                % 3 -> 6, no change
            end
        elseif nstate(i)==4
            if r1(i)<=g10*dt
                nstatetemp(i) = 5;
                % 4 -> 5, no change
            else
                xtemp(i) = x(i) + vr0*dt;
                % 4 -> 4, retrograde movement
                if (x(i)+leng(i)/2 - (floor(x(i)+leng(i)/2))) >= 0.5
                    k2 = (round(xtemp(i)+leng(i)/2)+1);
                else
                    k2 = (round(xtemp(i)+leng(i)/2));
                end             
                if nr~=k2   % moved by at least a bin, to the left
                    if nr-1<=axmin
                        nrp = nr+axmax;
                    else
                        nrp = nr;
                    end
                    
                    volume0(nl)      = volume0(nl) + 1;
                    volume0(nrp-1)   = volume0(nrp-1) - 1;
                    occupancy0(nl)   = occupancy0(nl) + 1;
                    occupancy0(nrp-1) = occupancy0(nrp-1) - 1;
                end
            end
        elseif nstate(i)==5
            if r1(i)<=g01*dt
                nstatetemp(i) = 4;
                xtemp(i) = x(i) + vr0*dt;
                % 5 -> 4, retrograde movement
                if (x(i)+leng(i)/2 - (floor(x(i)+leng(i)/2))) >= 0.5
                    k2 = (round(xtemp(i)+leng(i)/2)+1);
                else
                    k2 = (round(xtemp(i)+leng(i)/2));
                end   
                if nr~=k2   % moved by at least a bin, to the left
                    if nr-1<=axmin
                        nrp = nr+axmax;
                    else
                        nrp = nr;
                    end
                    
                    volume0(nl)       = volume0(nl) + 1;
                    volume0(nrp-1)    = volume0(nrp-1) - 1;
                    occupancy0(nl)    = occupancy0(nl) + 1;
                    occupancy0(nrp-1) = occupancy0(nrp-1) - 1;
                end
            elseif r1(i)<=(g01+goff)*dt
                nstatetemp(i) = 6;
                % 5 -> 6, lose one on track
                if nl+1>axmax
                    nlp = nl-axmax;
                else
                    nlp = nl;
                end
                if nr-1<=axmin
                    nrp = nr+axmax;
                else
                    nrp = nr;
                end
                
                for k=nlp+1:nrp-1 % Subtract from the old locations
                    if k>axmax
                        kk = k-axmax;
                    elseif k<1
                        kk = k+axmax;
                    else
                        kk=k;
                    end
                    occupancy0(kk) = occupancy0(kk)-1;
                end
            elseif r1(i)<=(g01+goff+gra)*dt
                nstatetemp(i) = 2;
                % 5 -> 2, no change
            end
        elseif nstate(i)==6
            if r1(i)<=gon*dt
                nstatetemp(i) = 5;
                % 6 -> 5, get one on track
                if nl+1>axmax
                    nlp = nl-axmax;
                else
                    nlp = nl;
                end
                
                if nr-1<=axmin
                    nrp = nr+axmax;
                else
                    nrp = nr;
                end
                
                for k=nlp+1:nrp-1 % Subtract from the old locations
                    if k>axmax
                        kk = k-axmax;
                    elseif k<1
                        kk = k+axmax;
                    else
                        kk=k;
                    end
                    occupancy0(kk) = occupancy0(kk)+1;
                end
            elseif r1(i)<=(gon+gra)*dt
                nstatetemp(i) = 3;
                % 6 -> 3, no change
            end
        end

        
        volume    = volume0;
        occupancy = occupancy0;
        
    end  % finish looping over NFs
    
    x = xtemp;
    nstate = nstatetemp;
    
    % enforce periodic boundary conditions
    for j = 1:nrunners
        if x(j)>=axmax
            xtemp(j) = x(j)-axmax;
        end
        
        if x(j)<=axmin
            xtemp(j) = x(j)+axmax;
        end
    end % finish second loop over NFs (for boundary conditions)
    x = xtemp;
    
    % Update total and on-track NF content at each location
    [volume]    = content(x,leng,axmax,nrunners); 
    [occupancy] = contenton(x,leng,axmax,nrunners,nstate); 
    
    % Update time
    temp   = temp+dt;
    
end

tnew = temp-dt;
xnew = x;
nstatenew    = nstate;
occupancynew = occupancy;
volumenew    = volume;



