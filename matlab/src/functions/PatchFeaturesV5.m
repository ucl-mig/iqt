function fts = PatchFeaturesV5(patch)

% 1:3 Eigenvalues of DT at patch centre
% 4:6 Principal eigenvector of DT at patch centre
% 7:10 Linearity, planarity, isotropy, trace.
% If n>0
% 11:13 Mean eigenvalues over central 3x3x3 patch
% 14:16 Mean principal direction over central 3x3x3 patch
% 17 Orientation dispersion over central 3x3x3 patch
% 18:21 Mean linearity, planarity, isotropy, trace over central 3x3x3 patch.
% If n>1
% 22:24 Mean eigenvalues over full patch
% 25:27 Mean principal direction over full patch
% 28 Orientation dispersion over full patch
% 29:32 Mean linearity, planarity, isotropy, trace over full patch.
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

ps=floor((length(patch)/6)^(1/3)+0.001);
rpatch = reshape(patch, [ps,ps,ps,6]);

esys = zeros(ps*ps*ps,10);
centre_ind = 0;
cpos = floor((ps+1)/2);
centre3x3_inds = [];
ind = 1;
for i=1:ps
    for j=1:ps
        for k=1:ps
            ldt = MakeDT_Matrix(rpatch(i,j,k,1),rpatch(i,j,k,2),rpatch(i,j,k,3),rpatch(i,j,k,4),rpatch(i,j,k,5),rpatch(i,j,k,6));
            [R, E] = eig(ldt);
            
            % Consistently oriented PD
            [av, bv] = eig(R(:,3)*R(:,3)');
            
            % Eigenvalues
            esys(ind,1:6) = [diag(E)', av(:,3)'];
            
            % Linearity
            esys(ind,7) = (E(3,3)-E(2,2))/E(3,3);
            % Planarity
            esys(ind,8) = (E(2,2)-E(1,1))/E(3,3);
            % Isotropy
            esys(ind,9) = E(1,1)/E(3,3);
            % Trace
            esys(ind,10) = E(1,1)+E(2,2)+E(3,3);

            if(i==cpos && j==cpos && k==cpos)
                centre_ind = ind;
            end
            if(abs(i-cpos)<=1 && abs(j-cpos)<=1 && abs(k-cpos)<=1)
                centre3x3_inds = [centre3x3_inds, ind];
            end
            ind = ind+1;
        end
    end
end

fts(1:10) = esys(centre_ind,:);

if(ps>1)
    mean_ev_3x3 = mean(esys(centre3x3_inds,1:3));

    pds_3x3 = esys(centre3x3_inds,4:6);
    Y = MeanDyadic(pds_3x3);
    [R E] = eig(Y);
    mean_pd_3x3 = R(:,3);
    dirconc_3x3 = E(3,3);

    mean_lpi_3x3 = mean(esys(centre3x3_inds,7:10));

    fts(11:13) = mean_ev_3x3;
    fts(14:16) = mean_pd_3x3;
    fts(17) = dirconc_3x3;
    fts(18:21) = mean_lpi_3x3;
end

if(ps>3)
    mean_ev_all = mean(esys(:,1:3));

    pds_all = esys(:,4:6);
    Y = MeanDyadic(pds_all);
    [R E] = eig(Y);
    mean_pd_all = R(:,3);
    dirconc_all = E(3,3);

    mean_fti_all = mean(esys(:,7:10));

    fts(22:24) = mean_ev_all;
    fts(25:27) = mean_pd_all;
    fts(28) = dirconc_all;
    fts(29:32) = mean_fti_all;

end

