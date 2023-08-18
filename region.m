function [r] = region(cov, time, tp_idx, cut_off)

    % Choose region
    if ( cov > cut_off)

        % R2 and R3 cov > 0.33
        if (time < tp_idx)
            r = 2;
        else
            r = 3;
        end
    else
        % R1 and R4 cov < 0.33
        if (time > tp_idx)
            r = 4;
        else
            r = 1;
        end
    end

end