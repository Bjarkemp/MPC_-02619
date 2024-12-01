function plots(t, u, Y)
    % Function to plot control inputs and tank heights over time.
    % Inputs:
    %   t - Time vector
    %   u - Control inputs matrix (size 2 x length(t))
    %   Y - Outputs matrix (size length(t) x 4)
    
    % Determine the final time from the time vector
    t_f = max(t);

    % Plot the control inputs F1 and F2 over time
    figure;
    for i = 1:2
        subplot(2, 1, i);
        plot(t, u(i, :), 'b', LineWidth = 2);
        xlabel('Time (s)');
        ylabel('Flow rate (cm^3/s)');
        % legend(['F' num2str(i)], 'Location', 'best');
        title(['Control Input F' num2str(i)]);
        xlim([0 t_f]);
    end

    % Plot the tank heights (outputs)
    figure;
    for i = 1:4
        subplot(2, 2, i);
        plot(t, Y(:, i), 'b', LineWidth = 2); % Assuming Y is given in row-wise time order
        xlabel('Time (s)');
        ylabel('Height (cm)');
        title(['Height in Tank ' num2str(i)]);
        % legend(['Tank ' num2str(i)], 'Location', 'best');
        xlim([0 t_f]);
    end
    hold on
end
