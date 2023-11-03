function sec = p2Sim
    
    load('Singapore.mat','data');
    lx        = length(data.B);
    data.int  = 5;
    data.tvec = [-75,linspace(data.Tres,365*3+1,data.int)];
    
    [data,dis,p2] = p2Params(data,'Covid Wildtype');
    
    xoptim     = [ones(lx,1), 1*ones(lx,1),0.5*ones(lx,1),1*ones(lx,1),1*ones(lx,1)]; %xoptim     = [ones(data.int*lx,1)];
    data.hw    = [zeros(data.int,lx)];
    %no testing in period 5 or after vaccine rollout (whatever comes first)
    data.imand = Inf;%no social distancing in period 1, but continues until end
    % data.imand is the index of the configuration with a stay at home
    % (SAH) order
  
    
    [data,f,g] = p2Run(data,dis,xoptim,p2);
    [cost,~]   = p2Cost(data,dis,p2,g);
    sec(1)     = sum(cost([3,6,7:10],:),'all');
    sec(2)     = sum(cost(3,:),'all');
    sec(3)     = sum(cost(6,:),'all');
    sec(4)     = sum(cost(7:10,:),'all');

    p2Plot(data,f,p2,g,cost);
    
end