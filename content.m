%% Count total number of NFs in space

function [volume] = content(x,leng,axmax,nrunners) 

volume = zeros(1,axmax);

for i = 1:nrunners
    
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

    % Count total NF content at each location
       for k = k1:k2
          
           if k>axmax
               kk = k-axmax;
           elseif k<1
               kk = k+axmax;
           else
               kk=k;
           end
           
           volume(kk) = volume(kk)+1;
       end

end
