function Pmean = mean_load_curve(P, delta_t)

intervals_per_day = 60*24/delta_t;

for i = 1:intervals_per_day
    
    ind = i:intervals_per_day:size(P,2)-intervals_per_day+i;
    
    Pmean(:,i) = mean(P(:,ind),2);
    
end

