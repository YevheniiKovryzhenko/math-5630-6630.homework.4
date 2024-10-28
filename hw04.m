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
            y = zeros(1,n_evals);
            y(:) = c(1);
            for i_pts = 1:n_evals
                for k = 2:n
                    cumprod = c(k);
                    for j = 1:k-1
                        cumprod = cumprod*(eval_(i_pts) - data(j,1));
                    end
                    y(i_pts) = y(i_pts)+cumprod;
                end
            end
        end


        function y = p2(data, eval)
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

            
            n = length(data);
            y = zeros(size(eval));

            % write your code here.
            
        end
    end
end
