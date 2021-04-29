% Developed by M Vito Hamza 07111640000154
% Based on Master Thesis by Annisaa Taradini
% DE Using YPEA107 DE
% function [batt]=inputbattDE
% clc;
% clear all;
% global data PL
tic %start timer
figure('Name','AI Worksheet');
CostFunction=@(P) ESinv(P);    % Cost Function

nVar=6;            % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin= [0 1];          % Lower Bound of Decision Variables
VarMax= [99999 30];          % Upper Bound of Decision Variables

% VarMin= data(:,4);          % Lower Bound of Decision Variables
% VarMax= data(:,5);          % Upper Bound of Decision Variables

%% mpoption OPF mba Dini
mpopt = mpoption('out.all',0, 'model', 'DC');

%% DE Parameters

MaxIt=3;      % Maximum Number of Iterations

nPop=4;        % Population Size

beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor

pCR=0.2;        % Crossover Probability

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
%     values = [1 2 3 4 5 6 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];
%     randidx = randperm(numel(values),2);
%     randvalues = values(randidx);
    pop(i).Position= [unifrnd(VarMin(1,1),VarMax(1,1),[1 4]) randperm(30,2)];  


%% Evalcost mba dini
    currentbatt=pop(i).Position;
%     case30mikrogridmaswince;
    rundcopf('case30mikrogridmaswince',mpopt,'','');
    pop(i).Cost=CostFunction(pop(i).Position);
   
  
    if pop(i).Cost<BestSol.Cost
    BestSol=pop(i);
    end 
    
end

BestCost=zeros(MaxIt,1);

for it=1:MaxIt
     for i=1:nPop
        subplot(2,2,[3,4]);
            xbus = [7 15 pop(i).Position(1,5) pop(i).Position(1,6)];
            xlim([1 30])
            title([it, i]);
            hold on
            scatter(xbus,pop(i).Position(1,1:4),'x');
            stem(xbus,BestSol.Position(1,1:4),'LineStyle','-.',...
                 'MarkerFaceColor','red',...
                 'MarkerEdgeColor','green');
            drawnow;
        
        x=pop(i).Position;
         
        A=randperm(nPop);
         
        A(A==i)=[];
        
        a=A(1)
        b=A(2)
        c=A(3)
        
        
%         beta=unifrnd(beta_min,beta_max);
        beta=unifrnd(beta_min,beta_max,VarSize);
%         ya=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        ya=pop(a).Position(1,1:4)+beta(1,1:4).*(pop(b).Position(1,1:4)-pop(c).Position(1,1:4));
        yb = round(pop(a).Position(1,5:6)+beta(1,5:6).*(pop(b).Position(1,5:6)-pop(c).Position(1,5:6)));
        ya = max(ya,VarMin(1,1)');
        ya = min(ya,VarMax(1,1)');
        yb = max(yb,VarMin(1,2)');
        yb = min(yb,VarMax(1,2)'); 
        round(yb);
%         yb(1,1) = 7;
%         yb(1,2) = 15;
%         yb(1,3) = 26;
        y = [ya yb];
        
         z=zeros(size(x));
         j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        cla(subplot(2,2,[3,4]));
%% Evalcost mba dini
        NewSol.Position=z;
        currentbatt=NewSol.Position;
%         case30mikrogridmaswince;
        rundcopf('case30mikrogridmaswince',mpopt,'','');
        NewSol.Cost=CostFunction(NewSol.Position);  
        
        if NewSol.Cost<pop(i).Cost 
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
               BestSol.Position
               BestSol.Cost
            end
        end
        
    end
      
    BestCost(it)=BestSol.Cost;

    for i = 1:it
        subplot(2,2,1);
            semilogy(BestCost, 'LineWidth', 2);
            xlabel('Iteration');
            ylabel('Best Cost');
            grid on;



        subplot(2,2,2);
            xbus = [7 15 BestSol.Position(1,5) BestSol.Position(1,6)];
            stem(xbus,BestSol.Position(1,1:4));
            xlim([0 30]);
            xlabel('Bus');
            ylabel('Capacity (Ah)');
            refreshdata;
            drawnow;
    end
end


%% Refresh Data to Best Solution
figure('Name','Best and Pop');
    xbus = [7 15 BestSol.Position(1,5) BestSol.Position(1,6)];
    hold on
    for sp1 = 1:nPop
    scatter(xbus,pop(sp1).Position(1,1:4),'x');
    drawnow
    end
    stem(xbus(1),BestSol.Position(1,1),'Color','g')
    stem(xbus(2),BestSol.Position(1,2),'Color','b')
    stem(xbus(3),BestSol.Position(1,3),'Color','r')
    stem(xbus(4),BestSol.Position(1,4),'Color','m')
    xlim([0 30]);
    xlabel('Bus');
    ylabel('Capacity (Ah)');

figure('Name','Semilogy');
semilogy(BestCost, 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Best Cost');
    grid on;
    
currentbatt = BestSol.Position
case30mikrogridmaswince;
rundcopf('case30mikrogridmaswince',mpopt,'','');

FINAL_BattInv = sum(0.209*72*BestSol.Position(1,1:4)*0.1295)       
FINAL_BestCostYear = FINAL_BattInv + Cost_Op*365

toc %stop timer