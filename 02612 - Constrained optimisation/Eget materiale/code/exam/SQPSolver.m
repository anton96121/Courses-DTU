function [x,z,Hist] = SQPSolver(x0,obj,con,l,u,cl,cu, options)
%{
This is the solver inteface for all the SQP algorithms. The interface has
 the input 'option' argument which functions as follows:
    Options:
    log - Do you want to log the process
        values: true/false

    method - Which SQP method
        values: {'SQP'        (Plain vanilla SQP)
                ,'SQP_ls'     (Line search SQP)
                ,'SQP_trust'  (Trust region SQP)}

    infesibility_handling - For infesibility handling SQP and SQP_ls 
        values: true/false

    precision - spans from 10^-1 to 10^-9
        values: Integer between 1 and 9

    subsolver - For SQP and SQP_ls without infesibility handling there is 
                a self made interior point algorithm. 
        values: 'own solver' or 'quadprog'

    trust_region - For SQP trust region one can set the inital trust region
        valuse: positive reals

    non_monotone - Activate or deactivatethe non monotone strategy for 
                   the line search SQP
        values: true/false

    penalty - Set the inital penalty
        values: reals above 100.


Created: 06.06.2021
Authors : Anton Ruby Larsen and Carl Frederik Grønvald
          IMM, Technical University of Denmark
%}
    if nargin <8
        options = struct();
    end
    if sum(strcmp(fieldnames(options), 'log')) == 1
        if islogical(options.log)
            log = options.log;
        else
            error('The option log should be a boolean')
        end
    else
        log = 0;
    end
    if sum(strcmp(fieldnames(options), 'method')) == 1
        if ismember(options.method, {'SQP', 'SQP_ls', 'SQP_trust'})
            method = options.method;
        else
            error('The option method should be one of {SQP, SQP_ls, SQP_trust}')
        end
    else
        method = 'SQP_ls';
    end
    if sum(strcmp(fieldnames(options), 'infesibility_handling')) == 1
       if islogical(options.infesibility_handling)
           infesibility_handling = options.infesibility_handling;
       else
           error('The options infesibility_handling should be a boolean')
       end
    else
        infesibility_handling = true;
    end
    if sum(strcmp(fieldnames(options), 'precision')) == 1
        if isnumeric(options.precision) && floor(options.precision)==options.precision && options.precision>0 && options.precision<6
            precision = options.precision;
        else
            error('The option precision should be a positive integer below 6.')
        end
    else
        precision = 3;
    end
    if sum(strcmp(fieldnames(options), 'subsolver')) == 1
        if ismember(options.subsolver , {'own solver', 'quadprog'})
            if strcmp(options.subsolver,'own solver')
                %true is own solver
                subsolver = true;
            else
                subsolver = false;
            end
        else
            msg = ['The option subsolver should be either be ',  char(39),  'own solver', char(39) ,' or ', char(39), 'quadprog' ,char(39)];
            error(msg)
        end
    else
        subsolver = true;
    end
    if sum(strcmp(fieldnames(options), 'trust_region')) == 1
        if isnumeric(options.trust_region) && options.trust_region>0
            trust_region = options.trust_region;
        else
            error('The option trust_region should be a positive real.')
        end
    else
        trust_region = 0.5;
    end
    if sum(strcmp(fieldnames(options), 'non_monotone')) == 1
        if islogical(options.non_monotone)
            non_monotone = options.non_monotone;
        else
            error('The option non_monotone should be a boolean.')
        end
    else
        non_monotone = true;
    end
        if sum(strcmp(fieldnames(options), 'penalty')) == 1
        if isnumeric(options.penalty) && options.penalty>100
            penalty = options.penalty;
        else
            error('The option penalty should be a real over 100.')
        end
    else
        penalty = 100;
    end
    
    
    
    if method == "SQP" && infesibility_handling==0
        try
        [x,z,Hist] = SQP(x0,obj,con,l,u,cl,cu,log,subsolver,precision);
        catch
            close all
            error('The program is infeasible. Try with infeasibility handling')
        end
        
    elseif method == "SQP_ls" && infesibility_handling==0
        try
        [x,z,Hist] = SQP_ls(x0,obj,con,l,u,cl,cu,log,subsolver,precision, non_monotone);
        catch
            close all
            error('The program is infeasible. Try with infeasibility handling')
        end
        
    elseif method == "SQP" && infesibility_handling==1
        [x,z,Hist] = SQP_infes(x0,obj,con,l,u,cl,cu,log,precision,penalty);
        
    elseif method == "SQP_ls" && infesibility_handling==1
        [x,z,Hist] = SQP_ls_infes(x0,obj,con,l,u,cl,cu,log,precision, non_monotone,penalty);
        
    elseif method == "SQP_trust"
        [x,z,Hist] = SQP_trust(x0,obj,con,l,u,cl,cu,log,precision, trust_region,penalty);
        
    end

end

