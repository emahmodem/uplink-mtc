clc;close all;

subject = 'Simulation Progress';

% scenarios = {  'CoverageThreshold_Pmin' , ...
%     'CoverageThreshold_SEPL' , ...
%     'CellDensity_Pmin'       , ...
%     'CellDensity_SEPL'       , ...
%     'M2MDensity_Pmin'       , ...
%     'M2MDensity_SEPL'       , ...
%     'ResourceBlocks_Pmin'       , ...
%     'PowerTruncationThreshold_Cell_Density'  ...
%     };

scenarios = {  
   
    'M2MDensity_SEPL'       
    
    
    };
for i = 1:length(scenarios)
    sim = scenarios{i};
    disp(char(datetime))
    try
        gen_coverage_results_mMTC_TWC(sim);
        send_report_via_email(subject,strcat(sim ,' has been finished at: ' , char(datetime)));
    catch err
        disp(err)
    end  
end
