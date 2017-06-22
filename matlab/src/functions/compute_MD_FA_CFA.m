function [md, fa, cfa ] = compute_MD_FA_CFA(dt)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

[XSIZE, YSIZE, ZSIZE, junk] = size(dt);
md = zeros(XSIZE,YSIZE,ZSIZE);
fa = zeros(XSIZE,YSIZE,ZSIZE);
cfa = zeros(XSIZE,YSIZE,ZSIZE,3);
parfor i=1:XSIZE
    for j=1:YSIZE
        for k=1:ZSIZE
            if(dt(i,j,k,3)>0)
                ldt = dt(i,j,k,:);
                md(i,j,k) = (ldt(3)+ldt(6) + ldt(8))/3;            
                ldt = MakeDT_Matrix(ldt(3),ldt(4),ldt(5),ldt(6),ldt(7),ldt(8));
                [R, E] = eig(ldt);
                fa(i,j,k) = sqrt(3*sum((diag(E) - mean(diag(E))).^2)/(2*sum(diag(E).^2)));
                cfa(i,j,k,:) = fa(i,j,k)*abs(R(:,3));
            end
        end
    end
end

md = squeeze(md);
fa = squeeze(fa);
cfa = squeeze(cfa);





