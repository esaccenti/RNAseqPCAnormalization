%% modify as needed
cd('..\RNAseqNORM\mfiles')
addpath('..\Addons\COVSCA\')
%%
clear all; close all; clc

%%
% Calculate agreement and conngrucne of score and loadings of different
% COVSCA models



FILES_IN{1} = "DATA_NEUROBLASTOMA_FULL_NORMALIZED_April2024.mat";
FILES_IN{2} = "DATA_TUMOR_FULL_NORMALIZED_April2024.mat";
FILES_IN{3} = "DATA_COLON_FULL_NORMALIZED_April2024.mat";

MaxR = 101; %Number of resamplings

for k = 1 : 3

    filein = fullfile('..\RNAseqNORM\data\data_full_matlab',FILES_IN{k});
    load(filein)

    [n,p] = size(DATA{1,1});

    for r = 1 : MaxR

        disp(sprintf('Repetition %i of %i', r, MaxR))

        idxr = randsample(p,500);

        for type = 1; % : 2

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


            RESULTS{k,r,1} = loadings;
            RESULTS{k,r,2} = scores;

        end

    end

end

%%
for k = 1 : 3

    jj = 0;

    % Set target

    TARGET = RESULTS{k,1,1};
    

    for i = 2 : MaxR
        for j = 2 : MaxR
            if i > j

                jj = jj + 1;

                Si = RESULTS{k,i,2};
                Sj = RESULTS{k,j,2};

                Si1 = Si(:,1);
                Si2 = Si(:,2);

                Sj1 = Sj(:,1);
                Sj2 = Sj(:,2);

                a = corr(Si1,Sj1);
                b = corr(Si1,Sj2);

                c = corr(Si2,Sj1);
                d = corr(Si2,Sj2);

                CS(k,jj) = max(abs([a b c d]));


                %%
                Si = RESULTS{k,i,1};
                Sj = RESULTS{k,j,1};

                Si1 = Si(:,1);
                Si2 = Si(:,2);

                Sj1 = Sj(:,1);
                Sj2 = Sj(:,2);

                a = corr(Si1,Sj1);
                b = corr(Si1,Sj2);

                c = corr(Si2,Sj1);
                d = corr(Si2,Sj2);

                CL(k,jj) = max(abs([a b c d]));

                %% first rotaion
                LtoRot1 = RESULTS{k,i,1};
                LtoRot2 = [LtoRot1(:,2),LtoRot1(:,2)];

                LtoRot3 = RESULTS{k,j,1};
                LtoRot4 = [LtoRot3(:,2),LtoRot3(:,2)];

                R1 = rotatefactors(LtoRot1,'Method','procrustes','Target',TARGET,'type','oblique');
                R2 = rotatefactors(LtoRot2,'Method','procrustes','Target',TARGET,'type','oblique');

                R3 = rotatefactors(LtoRot3,'Method','procrustes','Target',TARGET,'type','oblique');
                R4 = rotatefactors(LtoRot4,'Method','procrustes','Target',TARGET,'type','oblique');

                R = [R1 R2 R3 R4];

                [PF, comb] = TuckerPhiMat(R);


                CL(k,jj) = max(abs(PF));


            end
        end
    end

end

disp('Average Score correlation')
disp(mean(CS,2))


disp('Average Load congruence')
disp(mean(CL,2))


