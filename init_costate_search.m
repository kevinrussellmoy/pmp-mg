function [s_init_vec, SOC_end_vec] = init_costate_search(outage_len, pv_outage, ld_outage, ...
    x_max, x_min, inv_rating, inv_eff, ...
    socdot, P_batt_range_socdot, SOC_range_socdot, ...
    dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
    h, lambda_min, lambda_max, SOC_f_target)
    
deltaSOC_target = 0;
deltaSOC_tol = 0.01; % tolerance for accepting solution
MaxNumIter = 10;
numIter = 0;

s_min = lambda_min;
s_max = lambda_max;

    while numIter < MaxNumIter
    numIter = numIter+1;

        s_init = (s_min + s_max)/2;

        [~, SOC, ~] = mgpmpecm(outage_len, pv_outage, ld_outage, ...
                        x_max, x_min, inv_rating, inv_eff, ...
                        socdot, P_batt_range_socdot, SOC_range_socdot, ...
                        dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
                        h, s_init);

        SOC_end = SOC(end);
        deltaSOC = SOC_end - SOC_f_target;

        s_init_vec(numIter) = s_init;
        deltaSOC_vec(numIter) = deltaSOC;
        SOC_end_vec(numIter) = SOC_end;

%         disp(abs(deltaSOC - deltaSOC_target))
        disp(SOC_end)
        if abs(deltaSOC - deltaSOC_target) < deltaSOC_tol
            break  % exit "while" loop
        elseif deltaSOC > deltaSOC_target
            s_min = s_init;
        elseif deltaSOC < deltaSOC_target
            s_max = s_init;
        end

    end
    
    lambda_init = s_init_vec(end);
end