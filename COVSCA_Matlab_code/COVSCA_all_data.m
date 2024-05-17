%% Adapt as needed
cd('..\RNAseqNORM\mfiles')
addpath('..\Matlab\Addons\COVSCA\')

% Code for covsca available at
%https://github.com/esaccenti/covsca
%%
clear all; close all; clc

%%
dir_out_name = '..\RNAseqNORM\results\covsca';

FILES_OUT{1,1} = 'COVSCA_Blastoma_Full_Corr_April2024.xlsx';
FILES_OUT{1,2} = 'COVSCA_Blastoma_Full_Cov_April2024.xlsx';

FILES_OUT{2,1} = 'COVSCA_Tumor_Full_Corr_April2024.xlsx';
FILES_OUT{2,2} = 'COVSCA_Tumor_Full_Cov_April2024.xlsx';

FILES_OUT{3,1} = 'COVSCA_Colon_Full_Corr_April2024.xlsx';
FILES_OUT{3,2} = 'COVSCA_Colon_Full_Cov_April2024.xlsx';


FILES_IN{1} = "DATA_NEUROBLASTOMA_FULL_NORMALIZED_April2024.mat";
FILES_IN{2} = "DATA_TUMOR_FULL_NORMALIZED_April2024.mat";
FILES_IN{3} = "DATA_COLON_FULL_NORMALIZED_April2024.mat";

for k = 1 : 3

    filein = fullfile('..\RNAseqNORM\data\data_full_matlab',FILES_IN{k});
    load(filein)

    [n,p] = size(DATA{1,1});

    idxr = randsample(p,1000);

    for type = 1 : 2

        fileout = fullfile(dir_out_name,FILES_OUT{k,type});

        CORR = {};
        AllC = [];
        for f = 1:size(DATA,1)

            METHOD{f,1} =  DATA{f,2};

            X0 = DATA{f,1};

            X = X0(:,idxr);

            if type == 1 %corr

                C = corr(X);

            elseif type == 2 %cov

                C = cov(X);

            end

            CORR{f} = C;

            AllC = [AllC C];

            clear X

        end


        % Input parameters
        nanal = 10;  %Number of analysis

        %Rank and number of prototypes. 2 matrices of  this case, all of Rank 2
        %of rank 1
        Q= [1 1]';
        L = length(Q);


        % Run COVSCA
        [loadings, scores,fp,dys, func] = covsca(AllC,L,Q,9,1,nanal);

        % Fit percentages
        disp(fp)

        FP(k,type) = fp;

        %Correct GC name
        METHOD{4,1} = ' GCw';
        METHOD{5,1} = ' GCwb';

        % Save results
        sheet = 'covsca';
        HEAD = {'Method', 'CC1','CC2'};

        xlswrite(fileout,HEAD,sheet,'A1')
        xlswrite(fileout,METHOD,sheet,'A2')
        xlswrite(fileout,scores,sheet,'B2')




    end

end
