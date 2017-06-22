% Computes various features for each tensor in a tensor patch.  The
% features are:
%
% 1 Direction concentration of the tensor principal directions.
% 2-4 Average principle direction (x,y,z)
% 5 Mean first eigenvalue (L1).
% 6 Std L1.
% 7 Magnitude of L1 gradient over patch.
% 8 Absolute dot product of L1 gradient orientation with average 
% principle direction.
% Then the four features above for L1 repeated for each of the
% following:
%   9-12 Second eigenvalue L2.
%   13-16 Third eigenvalue L3.
%   17-20 Tensor linearity (L1-L2)/L1.
%   21-24 Tensor planarity (L2-L3)/L1.
%   25-28 Tensor isotropy L3/L1.
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%
function fts = PatchFeaturesV2(patch)

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
            % Must be a real rotation for quarternion representation so
            % reverse sense if necessary.  A rotation preserves cross
            % products so that R(v)xR(w) = R(vxw), whereas reflections
            % negate it.  This leads to the following test by setting v and
            % w to the x and y axes.
            if(abs(sum(cross(R(:,1), R(:,2)) - R(:,3)))>1E-6)
                R = -R;
            end
            q = qGetQ(R);
            [b, c, junk] = cart2sph(q(2), q(3), q(4));
            esys(ind,:) = [diag(E)', q(1), b, c];
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

mean_ev_3x3 = mean(esys(centre3x3_inds,1:3));
cov_ev_3x3 = cov(esys(centre3x3_inds,1:3));
mean_qu_3x3 = mean(esys(centre3x3_inds,4:6));
cov_qu_3x3 = cov(esys(centre3x3_inds,4:6));

fts(7:9) = mean_ev_3x3;
fts(10:12) = mean_qu_3x3;
fts(13:18) = cov_ev_3x3([1,2,3,5,6,9]);
fts(19:24) = cov_qu_3x3([1,2,3,5,6,9]);

if(ps>3)
    mean_ev_all = mean(esys(:,1:3));
    cov_ev_all = cov(esys(:,1:3));
    mean_qu_all = mean(esys(:,4:6));
    cov_qu_all = cov(esys(:,4:6));

    fts(25:27) = mean_ev_all;
    fts(28:30) = mean_qu_all;
    fts(31:36) = cov_ev_all([1,2,3,5,6,9]);
    fts(37:42) = cov_qu_all([1,2,3,5,6,9]);
end

