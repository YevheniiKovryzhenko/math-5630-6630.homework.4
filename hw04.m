% Author: Yevhenii Kovryzhenko / yzk0058@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw04


classdef hw04
    methods (Static)
        function y = p1(data, eval_)
            % Implement the divided difference method to interpolate the data points, then evaluate the polynomial at the given point.
            % :param data: a matrix of size n x 2, where the first column is the x values and the second column is the y values
            % :param eval: a column vector of x values at which to evaluate the polynomial
            % :return: a vector of y values corresponding to the evaluation points

            n = size(data, 1);
            
            D = zeros(n);
            D(:,1) = data(:,2);
            for k = 2:n
                for j = k:n
                    D(j,k) = (D(j,k-1)-D(j-1,k-1)) / (data(j, 1) - data(j-k+1, 1));
                end
            end
            
            c = diag(D);
            n_evals = length(eval_);
            y = zeros(1, n_evals);
            for k = 1:n_evals
                result = c(end);
                for i = n-1:-1:1
                    result = result * (eval_(k) - data(i,1)) + c(i);
                end
                y(k) = result;
            end
        end

    function y = p2(data, eval_)
        % for 6630 ONLY
        % Implement the divided difference method to interpolate the data 
        % points, then evaluate the polynomial at the given point.
        %
        % :param data: a cell array of size n x 1, each cell is a vector 
        %              like (x, y0, y1, y2, y3, ..., ym). 
        %              y0 = f(x), 
        %              y1 = f'(x), 
        %              y2 = f''(x),
        %              ... ,
        %              ym = f^{(m)}(x).
        % 
        %              Note, different cells may have different lengths.
        %
        % :param eval: a vector of x values at which to evaluate the polynomial
        % :return: a vector of y values corresponding to the evaluation points

        n = sum(cellfun(@(x) length(x) - 1, data));
        D = zeros(n);
        derivs = cell(1,n);
        
        x = zeros(1,n);
        idx = 0;
        for i = 1:length(data)
            idx = idx + 1;
            x(idx) = data{i}(1);
            D(idx,1) = data{i}(2);
            
            for ii = 3:length(data{i})
                idx = idx + 1;
                x(idx) = data{i}(1);
                D(idx,1) = data{i}(2);
                derivs{idx} = data{i}(3:end);
            end
        end

        
        for k = 2:n
            for j = k:n
                if x(j) == x(j-k+1)
                    D(j,k) = derivs{j}(k-1)/factorial(k-1);
                else
                    D(j,k) = (D(j,k-1)-D(j-1,k-1)) / (x(j) - x(j-k+1));
                end
            end
        end
    
        c = diag(D);    
        n_evals = length(eval_);
        y = zeros(1, n_evals);
        for k = 1:n_evals
            result = c(end);
            for i = n-1:-1:1
                result = result * (eval_(k) - x(i)) + c(i);
            end
            y(k) = result;
        end
    end



    end
end

