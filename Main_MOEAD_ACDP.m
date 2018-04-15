% MOEA/D-ACDP
function Main_MOEAD_ACDP()
clc;clear;tic;format compact;
problem_set    = {'LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12','LIRCMOP13','LIRCMOP14','CF1','CF2','CF3','CF4','CF5'}; % The number of test problems in this type
decision_dim   = 30 * ones(1,length(problem_set));

itrNum         = 30; % This program is run for 30 times


for p_id= 3:12
    
    itrCounter = 1;
    while itrCounter<=itrNum
        if(p_id>14)
            max_gen=1000;
        else
            max_gen    = 500;    % maximum generation number
        end
        popsize    = 300;    % population size
        
        %         T          = 0.1 * popsize; % the neighboor sizet
        T=20;
        n_r        = 2;
        n_dim      = decision_dim(p_id);     % the dimenstion of the decision vector
        delta      = 0.9; % the probability of selecting individuals from its neighboor.

        pcrossover = 1.0;    % crossover probability
        pmutation  = 1/n_dim;    % mutation probability
        de_factor  = 0.5;     % crossover distribution index
        mdi        = 20;     % mutation distribution index
        gen        = 1;      % the first generation
        lu         = decision_range(problem_set{p_id},n_dim)';
        fit_cv     = cmop_test(problem_set{p_id}); 
        rng(sum(100*clock));
%         pf=load(strcat(problem_set{p_id},'.pf'));
        % initialize the population
        x = ones(popsize, 1) * lu(1, :) + rand(popsize, n_dim) .* (ones(popsize, 1) * (lu(2, :) - lu(1, :))); %行为个体标识，列为变量标识
        % evaluate the population
        [f,cv] = fit_cv(x');
        Z           = min(f,[],2);
%         Z_feasible           = Inf * ones(size(f,1),1);  
        cv = overall_cv(cv);
        
        %ZQL initialize population and archive
        population=cell(popsize,1);
        M           = size(f,1); % the nubmer of objectives
        arch=[];
        for i=1:popsize
            if(cv(:,i)==0)
                population{i}.feasiblepop=[x(i,:),f(:,i)',cv(:,i)'];
                population{i}.infeasiblepop1=[];
                population{i}.infeasiblepop2=[];
            else
                population{i}.feasiblepop = [];
                population{i}.infeasiblepop1=[x(i,:),f(:,i)',cv(:,i)'];
                population{i}.infeasiblepop2=[x(i,:),f(:,i)',cv(:,i)'];
            end
            
            population{i}.dim=size([x(i,:),f(:,i)',cv(:,i)'],2);
            arch = archive(arch,[x(i,:),f(:,i)',cv(:,i)'],popsize,M);
        end
%         population=[x,f',cv'];
        %%%%%%%%%%%%% archive initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         fea=(population(:,end)==0);
%         arch=population(fea,:);

%         fea_ratio=sum(population(:,end)==0)/popsize;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (p_id<=12 || p_id>=15)
            lambda      =getlambda(popsize);% weight vectors
        else
            lambda      =load('W3D_300.dat');
        end
        B           = get_neighbors(popsize,lambda,T);
    
        while gen <= max_gen
%             tic;
            subproblem_selected  = randperm(popsize); 
            for i = 1:popsize
                % select and reproduce,P denotes the neighbors of chromosome i
                [Y,P] = select_and_reproduct(population,B,subproblem_selected(i),delta,n_dim,lu,pcrossover,de_factor,mdi,pmutation);
                [obj_Y ,con_Y] = fit_cv(Y');
                con_Y = overall_cv(con_Y);
                offspring = [Y,obj_Y',con_Y'];
                
                
                %update the value of Z
                
                
                Z = min(Z,obj_Y);
                
%                 if(con_Y==0)
%                     
%                 end
                %update the neighbors of chromosome i
                population = update_neighbors_NPOP(population,offspring,P,n_r,lambda,Z,n_dim,M);%% no fea_ratio and theta
            end
%             time(gen)=toc;
            %%%%%%%%%%%%%%%%% update archive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pop=[];
            for i=1:popsize
                pop=[pop;population{i}.feasiblepop];
            end
            if(~isempty(pop))
                arch = archive(arch,pop,popsize,M);
            end
            disp(sprintf('problem:%d, gen:%d',p_id,gen));
%             out=arch;
%             if isempty(out)
%                 out=[];
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             %  plot test
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             out_feasible=[];out_infeasible=[];
%             for ii=1:popsize
%                 if(~isempty([population{ii}.infeasiblepop1;population{ii}.infeasiblepop2]))
%                 out_infeasible=[out_infeasible;population{ii}.infeasiblepop1;population{ii}.infeasiblepop2];
%                 end
%                 if(~isempty(population{ii}.feasiblepop))
%                 out_feasible=[out_feasible;population{ii}.feasiblepop];
%                 end
%             end
%             if (p_id<=12 || p_id>=15)
%                         plot(pf(:,1),pf(:,2),'b+');
%                         if(~isempty(out))
%                         hold on;
%                         plot(out(:,n_dim+1),out(:,n_dim+2),'ro');
%                         end
%                         if(~isempty(out_infeasible))
%                             hold on;
%                             plot(out_infeasible(:,n_dim+1),out_infeasible(:,n_dim+2),'b*');
%                         end
%                         if(~isempty(out_feasible))
%                             hold on;
%                             plot(out_feasible(:,n_dim+1),out_feasible(:,n_dim+2),'g*');
%                         end
%                         hold off;
%                         xlabel('f1');ylabel('f2');
%                         title('External Archive');
%                         suptitle(strcat(problem_set{p_id},' , Generation: ',num2str(gen)));
%             else
%                         plot3(pf(:,1),pf(:,2),pf(:,3),'b+');
%                         hold on;
%                         plot3(out(:,n_dim+1),out(:,n_dim+2),out(:,n_dim+3),'ro');
%                         hold off;
%                         xlabel('f1');ylabel('f2');
%                         title('External Archive');
%                         suptitle(strcat(problem_set{p_id},' , Generation: ',num2str(gen)));
%             end
% 
%             drawnow;
            gen = gen + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % save the final results
        file_path = strcat(pwd,'/MOEAD_NPOP_Results/');
        if ~isdir(file_path)
            mkdir(file_path);
        end
        
        % show the percentage of the running
        
        out=arch;
        fit_obj = out(:,n_dim+1:end); x = out(:,1:n_dim);
        
        % objective and constraints fit_obj
        objName = strcat(file_path,problem_set{p_id},'_obj_cv_',num2str(itrCounter),'.dat');
        save(objName,'fit_obj','-ascii');
        decName = strcat(file_path,problem_set{p_id},'_dec_',num2str(itrCounter),'.dat');
        save(decName,'x','-ascii');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        itrCounter =itrCounter+1;
    end
end
% function y = teFit(fit,i)
%     lambda(i,:) = max(lambda(i,:),repmat(1e-6,size(lambda(i,:))));
%     y =max(1./lambda(i,:).*abs(fit - z')); 
% end
end


