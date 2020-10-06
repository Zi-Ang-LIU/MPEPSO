%***************************************************************************************************************************************
%  MPEPSO 
%  Version: MPEPSO 1.0
%  Programmed by Ziang Liu at Osaka University 2020
%  Email: ziang@inulab.sys.es.osaka-u.ac.jp;
%  Date: 2020/10/5

%  Paper: 
%  [1] Ziang Liu and Tatsushi Nishi, "Multi-population Ensemble Particle Swarm Optimizer for Engineering Design Problems"
%      Mathematical problems in engineering, 2020. https://doi.org/10.1155/2020/1450985

%  References: 
%  [1] Lynn N, Suganthan P N. Ensemble particle swarm optimizer[J]. Applied Soft Computing, 2017, 55: 533-548.
%  [2] Wu G, Mallipeddi R, Suganthan P N, et al. Differential evolution with multi-population based ensemble of mutation strategies[J]. 
%      Information Sciences, 2016, 329: 329-345.
%  [3] Wu G, Shen X, Li H, et al. Ensemble of differential evolution variants[J]. Information Sciences, 2018, 423: 172-186.
%***************************************************************************************************************************************

%% Clear workspace and command window
clc;
clear;

rand('state', sum(100*clock));

%mex cec14_func.cpp -DWINDOWS
% f = cec14_func(x,func_num); here x is a D*pop_size matrix.

%D=10, 30, 50, 100; Runs / problem: 51;
%MaxFES: 10000*D (Max_FES for 10D = 100000; for 30D = 300000; for 50D = 500000; for 100D = 1000000)


%% Initialize parameters

% run 30 times
run=30;

% dimension number
dimension=10;

% search range
range=[-100 100];

% optimal solution
optima = cumsum(100*ones(1,30));

% Preallocating a matrix for solutions
solution=zeros(30,run);



% parameter setting
indexBestLN = 1;
popsize=50;
leastSelectionPro=0.1;
genForChange=10;

% number of function evaluations
max_FES=10000*dimension;
% number of iterations
max_iteration=max_FES/popsize;  
% Preallocating a matrix for convergence graphs
con_graph=ones(run,max_FES,30);

%% Iteration

for func_num=1:30
    
    for i=1:run
        
        %MPEPSO
        [position,value,convergence] = MPEPSO(popsize,range,dimension,max_iteration,max_FES,func_num,genForChange,indexBestLN,leastSelectionPro);

        % data for convergence graph
        con_graph(i,:,func_num)=convergence(1:max_FES);
        
        % record error
        solution(func_num,i) = value-optima(func_num);
        
    end
    
    % mean error
    m = mean(solution(func_num,:),2);
    
    % standard deviation
    s = std(solution(func_num,:),0,2);
    
    % output the results
    fprintf('Func_%d\n Mean:\t%d\n Std:\t%d\n', func_num, m, s);
    
end

%% Save

% algo_name='MPEPSO_';
% file_name= [algo_name,num2str(popsize),'pop_',num2str(dimension),'D_',num2str(genForChange),'LP_',num2str(indexBestLN),'indexBestLN_',num2str(leastSelectionPro),'leastSelectionPro_',num2str(run),'run','.mat'];
% save(file_name,'solution')
% 
% 
% convergence_name= [algo_name,'convergence_',num2str(popsize),'pop_',num2str(dimension),'D_',num2str(genForChange),'LP_',num2str(indexBestLN),'indexBestLN_',num2str(leastSelectionPro),'leastSelectionPro_',num2str(run),'run','.mat'];
% save(convergence_name,'con_graph', '-v7.3')

%% convergence graphs
% for i=1:30
%     figure
%     semilogy(median(con_graph(:,:,i)));
% end
