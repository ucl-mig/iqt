function fts = PatchFeaturesV4(patch)

% 1:3 Eigenvalues of DT at patch centre
% 4:6 Principal eigenvector of DT at patch centre
% If n>0
% 7:9 Mean eigenvalues over central 3x3x3 patch
% 10:12 Mean principal direction over central 3x3x3 patch
% 13 Orientation dispersion over central 3x3x3 patch
% If n>1
% 14:16 Mean eigenvalues over full patch
% 17:19 Mean principal direction over full patch
% 20 Orientation dispersion over full patch

ps=floor((length(patch)/6)^(1/3)+0.001);
rpatch = reshape(patch, [ps,ps,ps,6]);

esys = zeros(ps*ps*ps,6);
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
            
            esys(ind,:) = [diag(E)', av(:,3)'];
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

fts(1:6) = esys(centre_ind,:);

if(ps>1)
    mean_ev_3x3 = mean(esys(centre3x3_inds,1:3));

    pds_3x3 = esys(centre3x3_inds,4:6);
    Y = MeanDyadic(pds_3x3);
    [R E] = eig(Y);
    mean_pd_3x3 = R(:,3);
    dirconc_3x3 = E(3,3);

    fts(7:9) = mean_ev_3x3;
    fts(10:12) = mean_pd_3x3;
    fts(13) = dirconc_3x3;
end

if(ps>3)
    mean_ev_all = mean(esys(:,1:3));

    pds_all = esys(:,4:6);
    Y = MeanDyadic(pds_all);
    [R E] = eig(Y);
    mean_pd_all = R(:,3);
    dirconc_all = E(3,3);

    fts(14:16) = mean_ev_all;
    fts(17:19) = mean_pd_all;
    fts(20) = dirconc_all;

end

