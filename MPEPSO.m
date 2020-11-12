function [position,value,convergence] = MPEPSO(popsize,range,dimension,max_iteration,max_FES,func_num,genForChange,indexBestLN,leastSelectionPro)
%***************************************************************************************************************************************
%  MPEPSO
%  Version: MPEPSO 1.0
%  Programmed by Ziang Liu at Osaka University 2020
%  Email: ziang@inulab.sys.es.osaka-u.ac.jp;
%  Date: 2020/10/5

%  Paper:
%  [1] Ziang Liu and Tatsushi Nishi, "Multi-population Ensemble Particle Swarm Optimizer for Engineering Design Problems"
%      Mathematical Problems in Engineering, 2020. https://doi.org/10.1155/2020/1450985

%  References:
%  [1] Lynn N, Suganthan P N. Ensemble particle swarm optimizer[J]. Applied Soft Computing, 2017, 55: 533-548.
%  [2] Wu G, Mallipeddi R, Suganthan P N, et al. Differential evolution with multi-population based ensemble of mutation strategies[J].
%      Information Sciences, 2016, 329: 329-345.
%  [3] Wu G, Shen X, Li H, et al. Ensemble of differential evolution variants[J]. Information Sciences, 2018, 423: 172-186.
%***************************************************************************************************************************************

%% initialization
FES = 0;
strategy_num=3;
arrayGbestChange = ones(1,strategy_num);
arrayGbestChangeRate = zeros(1,strategy_num);
%allocation count
numViaLN = zeros(1,strategy_num);

if size(range,1)==1
    range_min=range(1)*ones(popsize,dimension);
    range_max=range(2)*ones(popsize,dimension);
else
    range_min=range(1,:);
    range_min=repmat(range_min,popsize,1);
    range_max=range(2,:);
    range_max=repmat(range_max,popsize,1);
end

interval = range_max-range_min;
v_max=interval*0.5;
v_min=-v_max;

% Initialize positions and velocities
pos = range_min+ interval.*rand(popsize,dimension);
vel =v_min+(v_max-v_min).*rand(popsize,dimension);

% function evaluation (CEC 2014)
val = (cec14_func(pos',func_num))';

% update the number of function evaluations
FES=FES+popsize;

% calculate gbest
[gbest_val,g_index]=min(val);
gbest_pos=pos(g_index,:);

convergence(1:FES)=gbest_val;

% calculate pbest
pbest_pos=pos;
pbest_val=val';

permutation = randperm(popsize);
arrayThird= permutation(1:leastSelectionPro*popsize);
arraySecond = permutation(leastSelectionPro*popsize+1: 2*leastSelectionPro*popsize);
arrayFirst = permutation(2*leastSelectionPro*popsize+1:end);

%% Initialize LDWPSO
% LDWPSO
w2=0.9-(1:max_iteration)*(0.7/max_iteration);
c2_1=2.5-(1:max_iteration)*2/max_iteration;
c2_2=0.5+(1:max_iteration)*2/max_iteration;

%% Initialize UPSO
% UPSO
c6_1=2.5-(1:max_iteration)*2/max_iteration;   %acceleration constants
c6_2=0.5+(1:max_iteration)*2/max_iteration;
w_6=0.9-(1:max_iteration)*(0.7/max_iteration);
u_6=0.1;
mu_6=0;
sigma_6=0.01;
nor=normrnd(mu_6,sigma_6,max_FES,dimension);

neighbor_6(1,:)=[popsize,2];
for i=2:(popsize-1)
    neighbor_6(i,:)=[i-1,i+1];
end
neighbor_6(popsize,:)=[popsize-1,1];

%% Initialize CLPSO
m=5;
obj_func_slope=zeros(popsize,1);
fri_best_pos=zeros(popsize,dimension);
fri_best=(1:popsize)'*ones(1,dimension);
j=0:(1/(popsize-1)):1;
j=j*10;
Pc=ones(dimension,1)*(0.0+((0.5).*(exp(j)-exp(j(1)))./(exp(j(popsize))-exp(j(1)))));

for i=1:popsize
    fri_best(i,:)=i*ones(1,dimension);
    friend1=ceil(popsize*rand(1,dimension));
    friend2=ceil(popsize*rand(1,dimension));
    friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
    toss=ceil(rand(1,dimension)-Pc(:,i)');
    if toss==ones(1,dimension)
        temp_index=randperm(dimension);
        toss(1,temp_index(1))=0;
        clear temp_index;
    end
    fri_best(i,:)=(1-toss).*friend+toss.*fri_best(i,:);
    for d=1:dimension
        fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
    end
end

%% Iteration

% number of iterations
k=0;

% allocation rate
rateViaLN=zeros(max_iteration, 3);

while k<=max_iteration && FES<=max_FES
    
    k=k+1;
    gbest_pos_temp=repmat(gbest_pos,popsize,1);
    
    if mod(k,genForChange) == 0
        arrayGbestChangeRate(1) = arrayGbestChange(1)/length(arrayFirst);
        arrayGbestChangeRate(2) = arrayGbestChange(2)/length(arraySecond);
        arrayGbestChangeRate(3) = arrayGbestChange(3)/length(arrayThird);
        [~,indexBestLN]=max(arrayGbestChangeRate);
        arrayGbestChange = [1,1,1];
        arrayGbestChangeRate =  [0,0,0];
    end
    
    permutation = randperm(popsize);
    if indexBestLN == 1
        arrayThird= permutation(1:leastSelectionPro*popsize);
        arraySecond = permutation(leastSelectionPro*popsize+1: 2*leastSelectionPro*popsize);
        arrayFirst = permutation(2*leastSelectionPro*popsize+1:end);
        numViaLN(1) = numViaLN(1) + 1;
    elseif indexBestLN == 2
        arrayThird = permutation(1:leastSelectionPro*popsize);
        arrayFirst = permutation(leastSelectionPro*popsize+1: 2*leastSelectionPro*popsize);
        arraySecond  = permutation(2*leastSelectionPro*popsize+1:end);
        numViaLN(2) = numViaLN(2) + 1;
    elseif indexBestLN == 3
        arrayFirst = permutation(1:leastSelectionPro*popsize);
        arraySecond = permutation(leastSelectionPro*popsize+1: 2*leastSelectionPro*popsize);
        arrayThird  = permutation(2*leastSelectionPro*popsize+1:end);
        numViaLN(3) = numViaLN(3) + 1;
    end
    rateViaLN(k,:) = numViaLN/sum(numViaLN);%allocation rate
    
    % LDWPSO
    popsize1= length(arrayFirst);
    vel(arrayFirst,:)=w2(k).*vel(arrayFirst,:)+(c2_1(k).*rand(popsize1,dimension).*(pbest_pos(arrayFirst,:)-pos(arrayFirst,:))+c2_2(k).*rand(popsize1,dimension).*(gbest_pos_temp(arrayFirst,:)-pos(arrayFirst,:)));
    
    % UPSO
    for i=arraySecond
        [~,tmpid]=min(pbest_val(neighbor_6(i,:)));
        % aa1(i,:)=c6_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+c6_2(k).*rand(1,dimension).*(gbest_pos-pos(i,:));
        vel1(i,:)=w_6(k).*vel(i,:)+( c6_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+c6_2(k).*rand(1,dimension).*(gbest_pos-pos(i,:)) );
        % aa2(i,:)=c6_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+c6_2(k).*rand(1,dimension).*(pbest_pos(neighbor_6(i,tmpid),:)-pos(i,:));
        vel2(i,:)=w_6(k).*vel(i,:)+( c6_1(k).*rand(1,dimension).*(pbest_pos(i,:)-pos(i,:))+c6_2(k).*rand(1,dimension).*(pbest_pos(neighbor_6(i,tmpid),:)-pos(i,:)) );
        r3_6=nor(FES,:);
        vel(i,:)=r3_6.*u_6.*vel1(i,:)+(1-u_6).*vel2(i,:);
    end
    
    % CLPSO
    for i=arrayThird
        delta(i,:)=(c2_1(k).*rand(1,dimension).*(fri_best_pos(i,:)-pos(i,:)));
        vel(i,:)=w2(k)*vel(i,:)+delta(i,:);
        if obj_func_slope(i)>m
            fri_best(i,:)=i*ones(1,dimension);
            friend1=(ceil(popsize*rand(1,dimension)));
            friend2=(ceil(popsize*rand(1,dimension)));
            friend=(pbest_val(friend1)<pbest_val(friend2)).*friend1+(pbest_val(friend1)>=pbest_val(friend2)).*friend2;
            toss=ceil(rand(1,dimension)-Pc(:,i)');
            if toss==ones(1,dimension)
                temp_index=randperm(dimension);
                toss(1,temp_index(1))=0;
                clear temp_index;
            end
            fri_best(i,:)=(1-toss).*friend+toss.*fri_best(i,:);
            for d=1:dimension
                fri_best_pos(i,d)=pbest_pos(fri_best(i,d),d);
            end
            obj_func_slope(i)=0;
        end
    end
    
    %% Update positions
    vel=((vel<v_min).*v_min)+((vel>v_max).*v_max)+(((vel<v_max)&(vel>v_min)).*vel);
    pos=pos+vel;
    
    for i=1:popsize
        if (sum(pos(i,:)>range_max(i,:))+sum(pos(i,:)<range_min(i,:))==0)
            val(i) = (cec14_func(pos(i,:)',func_num))';
            FES=FES+1;
            
            convergence(FES)=gbest_val;
            
            if FES>=max_FES
                break;
            end
            
            % update pbest
            if  val(i)<pbest_val(i)
                if ismember(i,arrayFirst)
                    arrayGbestChange(1) = arrayGbestChange(1) + pbest_val(i)-val(i);
                elseif ismember(i,arraySecond)
                    arrayGbestChange(2) = arrayGbestChange(2) + pbest_val(i)-val(i);
                elseif ismember(i,arrayThird)
                    arrayGbestChange(3) = arrayGbestChange(3) + pbest_val(i)-val(i);
                end
                pbest_pos(i,:)=pos(i,:);
                pbest_val(i)=val(i);
                obj_func_slope(i)=0;
            else
                obj_func_slope(i)=obj_func_slope(i)+1;
            end
            
            % update gbest
            if  pbest_val(i)<gbest_val
                gbest_pos=pbest_pos(i,:);
                gbest_val=pbest_val(i);
            end
        end
    end
    
    if FES>=max_FES
        break;
    end
    
    if (k==max_iteration)&&(FES<max_FES)
        k=k-1;
    end
    
end

position=gbest_pos;
value=gbest_val;
% plot(rateViaLN)

end
