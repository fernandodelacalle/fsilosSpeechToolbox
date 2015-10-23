function [SE] = se3dNew(t_fin, t_ini, time_step, frequencySize)
    
    [M_temp , t] = temporalMaskingNew(t_fin,t_ini, time_step);
    [M_frec, ~] = frecuencyMaskingNew(frequencySize, 0);

    SE = zeros(size(t,2), size(M_frec,2));
    for j = 1:size(SE,2)
        SE(:,j) = M_temp;
    end

    for i = 1:size(SE,1)   
        SE(i,:) = (SE(i,:) .* M_frec);
    end
end

