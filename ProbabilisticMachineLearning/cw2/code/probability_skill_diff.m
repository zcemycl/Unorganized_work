function [p1,p2,p3] = probability_skill_diff(m1,m2,o1,o2, mode)
    if mode == 1
        m_1 = m2;
        m_2 = m1;
        o_1 = o2;
        o_2 = o1;
    elseif mode == 0 
        m_1 = m1;
        m_2 = m2;
        o_1 = o1;
        o_2 = o2;
    end
    p1 = normcdf((m_2-m_1)/sqrt(o1^2+o2^2),0,1);
    p2 = (normcdf((m_2-m_1)/sqrt(1+o1^2+o2^2),0,1));
%     p3 = normcdf(0,0,1)*(normcdf((m2-m1)/o1,0,1));
    
end