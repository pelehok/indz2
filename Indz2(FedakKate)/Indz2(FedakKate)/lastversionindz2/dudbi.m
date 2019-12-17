function dudb = dudbi(t,params)
if params.n==0 
    count=1;
else
    step=params.step;
    count=floor(t/step)+1;
    count=min(count,params.n);
end
t1=params.t0+(count-1)*step;
t2=t1+step;
if (params.n==0) 
    dudb=1; 
else
if (count==params.k) 
    dudb = (t2-t)/step;       
elseif (count+1==params.k && count ~= t/step+1)
    dudb = (t-t1)/step;         
else
    dudb=0;
end
end
end

