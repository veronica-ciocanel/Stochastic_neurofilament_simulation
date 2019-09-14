%% Count number of on-track NFs in space

function [occupancy] = contenton(x,leng,axmax,nrunners,nstate)

occupancy = zeros(1,axmax);

for i = 1:nrunners
    
    if nstate(i)==1 || nstate(i)==2 || nstate(i)==4 || nstate(i)==5
        
        nr = round(x(i)+leng(i)/2); % bin of the right end of nf_i
        nl = round(x(i)-leng(i)/2); % bin of the left end of nf_i
        
        % Round bins covered by the NF based on how far it extends into bins  
        if (x(i)+leng(i)/2 - (floor(x(i)+leng(i)/2))) >= 0.5
            k1 = nl + 1;
            k2 = nr;
        else
            k1 = nl;
            k2 = nr - 1;
        end
        
        % Count on-track NF content at each location
        for k = k1:k2
            
            if k>axmax
                kk = k-axmax;
            elseif k<1
                kk = k+axmax;
            else
                kk=k;
            end
            
            occupancy(kk) = occupancy(kk)+1;
        end
    end  
    
end

