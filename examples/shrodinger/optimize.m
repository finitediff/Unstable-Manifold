function [outFunction] = optimize(inFunction)
%OPTIMIZE Optimizes a function passed in by saving parameters.
%   It takes in a function and returns a function.  When the returned
%   function is called, it records the parameters and saves the output to a
%   file.  Thus, optimized functions will perform faster by loading
%   previously saved values.

    outFunction = @ (varargin) optimizeHelper(inFunction, varargin);
end


function solution = optimizeHelper(varargin)
    inFunction = varargin{1};
    parameters = varargin{2};
    parameterString = "";
    
    % Get the file name.
    for i=1:size(parameters,2)
        temp = parameters{i};
        if isa(temp, 'string')
                
        parameterString = parameterString + temp;
        else
            if isa(temp, 'double') && size(temp,1) == 1 && size(temp,2) == 1
                parameterString = parameterString + num2str(temp, "%.2f");

            end
        end
 
        parameterString = parameterString + "__";
    end
    
    fileStr = strcat(pwd,"/generated/", func2str(inFunction)) +"__"+ ...
        parameterString;
    fileStr = fileStr.replace(".","p") + ".mat";
    
    % If the file already exists
    if exist(fileStr, 'file') == 2
        fprintf(" Loading file");
        ld = load(fileStr);
        solution = ld.solution;
        
    % If it doesn't exist
    else
        fprintf(" Generating and Saving");
        solution = inFunction(parameters{:});
        save(fileStr, 'solution');
    end
    
    
end