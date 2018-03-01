%% Question 1
% Define interval [xl xu]
xl = 0.4; xu = 0.48;
% Initialize midpoint
f_xr = 1;
% Initialize counter
i = 1;

% While the function is above tolerance
while abs(f_xr) >= 1e-5
    % Find midpoint of interval
    xr = (xl + xu)/2;
    % Evaluate function at xr, xl, and xu
    f_xr = tan(pi*xr) - xr - 6;
    f_xl = tan(pi*xl) - xl - 6;
    f_xu = tan(pi*xu) - xu - 6;
    % Relative error from previous step
    error = xr - 0.4510472588;
    % Store values for post-processing
    tab(i,:) = [i xl xu xr f_xr error];
    % Find where the function changes sign and assign new interval
    if f_xu > 0 && f_xr < 0     % if sign changes in this domain
        xl = xr;                % Assign xl as a bound
    elseif f_xr > 0 && f_xl < 0 % if sign changes in this domain
        xu = xr;                % Assign xu as a bound
    end
    i = i + 1;  % Update counter
end

% Plotting
x = linspace(0,5,1000);
plot(x,(tan(pi*x) - x - 6));
axis([0 5 -10 10]); title('Function Plot'); ylabel('y'); xlabel('x');


%% Question 2
%% Part a.)
% Define initial guess
xi_1 = 0.48;
% Initialize f(xi)
f_x = 1;
% Initialize counter
k = 1;

% While the function is above tolerance
while abs(f_x) >= 1e-6
    % Limit to 25 iterations
    if k > 25
        return
    end
    % Report table values
    tab2(k,1:2) = [k xi_1];
    % Assign new value x_i+1 to xi
    xi = xi_1;
    % Find value of f(xi)
    f_x = tan(pi*xi) - xi - 6; 
    % Find value of df/dx at xi
    df_dx = pi*sec(pi*xi)^2 - 1;
    % Use N-R equation
    xi_1 = xi - f_x/df_dx;
    % Report table values
    tab2(k,3:4) = [xi_1 xi_1-0.4510472588];
    % Update counter
    k = k + 1;
end


%% Part b
% Define initial guesses
x0 = 0.54;
x1 = 0.48;
% Initialize counter and f(x+1)
k = 1;
f_x1 = 1;
while abs(f_x1) >= 1e-5
    % Limit to 25 iterations
    if k > 25
        return
    end
    % Find value of f(x1)
    f_x1 = tan(pi*x1) - x1 - 6; 
    % Find value of f(x0)
    f_x0 = tan(pi*x0) - x0 - 6;
    % Use Secant equation
    x2 = x1 - f_x1 * (x1-x0)/(f_x1-f_x0);
    tab3(k,:) = [k x1 x2 x2-0.4510472588];
    % Update variables for next iteration
    x0 = x1;
    x1 = x2;
    % Update counter
    k = k + 1;
end

% Plotting for part e.)
figure(1)
semilogy(tab2(:,1)',tab2(:,4)')
hold on
semilogy(tab3(:,1)',tab3(:,4)')
title('Relative Error of Solution'); xlabel('Iteration');
ylabel('Relative Error');
legend('Newton-Raphson','Secant');

%% Question 3
% Define interval [xl xu]
xl = 0.4; xu = 0.48;
% Initialize midpoint
f_xr = 1;
% Initialize counter
i = 1;

% While the function is above tolerance
while abs(f_xr) >= 1e-5
    % Store old f_xr
    f_xr_old = f_xr;
    % Find midpoint of interval
    xr = (xl + xu)/2;
    % Evaluate function at xr, xl, and xu
    f_xr = tan(pi*xr) - xr - 6;
    f_xl = tan(pi*xl) - xl - 6;
    f_xu = tan(pi*xu) - xu - 6;
    % Relative error from previous step
    error = f_xr_old - f_xr;
    % Store values for post-processing
    tab(i,:) = [i xl xu xr f_xr error];
    % Find where the function changes sign and assign new interval
    if f_xu > 0 && f_xr < 0
        xl = xr;
    elseif f_xr > 0 && f_xl < 0
        xu = xr;
    end
    i = i + 1;
end

