function [association loads] = associate_nodes(locations,mode)
switch(mode)
    case 'AGGREGATION'
        distances_H2H_to_BS = pdist2(locations.BS,locations.H2H,'euclidean') ;
        [~ , server_h2h] = min(distances_H2H_to_BS);
        association.H2H_to_BS = server_h2h;
        
        distances_CH_to_BS = pdist2(locations.BS,locations.CH,'euclidean') ;
        [~ , server_ch] = min(distances_CH_to_BS);
        association.CH_to_BS = server_ch;
        
        aggregators = [locations.H2H ; locations.CH];
        N_H2H = size(locations.H2H ,1);
        N_CH = size(locations.CH ,1);
        distances_M2M_to_H2H = pdist2(aggregators,locations.M2M,'euclidean') ;
        [~ , server_m2m] = min(distances_M2M_to_H2H);
        
        Im2m_h2h = find(server_m2m <= N_H2H);
        Im2m_ch = find(server_m2m > N_H2H);
        Ih2h = server_m2m(Im2m_h2h);
        Ich = server_m2m(Im2m_ch) - N_H2H  ;
        
        association.M2M_to_H2H = [Im2m_h2h' Ih2h'];
        association.M2M_to_CH = [Im2m_ch' Ich'];
        
        for i = 1:N_H2H
            loads.H2H(i) = sum(association.M2M_to_H2H(:,2) == i);
        end
        
        for i = 1:N_CH
            loads.CH(i) = sum(association.M2M_to_CH(:,2) == i);
        end
        
        for i = 1:size(locations.BS ,1)
            loads.UE_H2H(i) = sum(association.H2H_to_BS == i);
            loads.UE_CH(i) = sum(association.CH_to_BS == i);
        end
        
    case 'DIRECT'
        
        distances_H2H_to_BS = pdist2(locations.BS,locations.H2H,'euclidean') ;
        [~ , server_h2h] = min(distances_H2H_to_BS);
        association.H2H_to_BS = server_h2h;
        
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [~ , server_m2m] = min(distances_M2M_to_BS);
        association.M2M_to_BS = server_m2m;
        
        N_BS = size(locations.BS ,1);
        for i = 1:N_BS
            loads.BS(i) = sum(association.M2M_to_BS == i);
        end
end
end