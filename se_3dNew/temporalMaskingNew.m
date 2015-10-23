function [ M, t ] = temporalMaskingNew(t_fin,t_ini, step)

    maxdB=60;
    t_pos = 0:step:t_fin;
    M_pos  = ellipse(t_fin,maxdB,t_pos);
    m_pos = 1/t_fin;
    M_pos = (-m_pos*t_pos+1)*maxdB;

    t_pre = t_ini:step:0;

    M_pre = ellipse(t_ini,maxdB,t_pre);
    
    t = [t_pre t_pos(2:end)];
    M = [M_pre M_pos(2:end)];
end