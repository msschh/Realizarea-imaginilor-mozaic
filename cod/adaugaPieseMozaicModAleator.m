function imgMozaic = adaugaPieseMozaicModAleator( params )
imgMozaic = uint8(zeros(size(params.imgReferintaRedimensionata)));
[H,W,~,N] = size(params.pieseMozaic);
[h,w,c] = size(params.imgReferintaRedimensionata);
p = h - H + 1;
q = w - W + 1;
a = 1:p;
b = 1:q;
[X,Y] = meshgrid(a, b);
pairs = [X(:) Y(:)];
perm = randperm(p * q);

if c == 3
    medii = zeros(N, 3);
    for i = 1:N
        medii(i, 1) = mean2(params.pieseMozaic(:, :, 1, i));
        medii(i, 2) = mean2(params.pieseMozaic(:, :, 2, i));
        medii(i, 3) = mean2(params.pieseMozaic(:, :, 3, i));
    end
    for l = 1:(p * q)
        i = pairs(perm(l), 1);
        j = pairs(perm(l), 2);
        if (imgMozaic(i, j, 1) == 0 && imgMozaic(i, j, 2) == 0 && imgMozaic(i, j, 3) == 0) || i == p || j == q 
            mR = mean2(params.imgReferintaRedimensionata(i : i + H - 1,j : j + W - 1, 1));
            mG = mean2(params.imgReferintaRedimensionata(i : i + H - 1,j : j + W - 1, 2));
            mB = mean2(params.imgReferintaRedimensionata(i : i + H - 1,j : j + W - 1, 3));
            distances = zeros(N);
            for k = 1:N
                distances(k) = (mR - medii(k,1)).^2 + (mG - medii(k,2)).^2 + (mB - medii(k,3)).^2;
            end
            [~ , indiceAles] = min(distances(:,1));
            imgMozaic(i : i + H - 1,j : j + W - 1, :) = params.pieseMozaic(:,:,:,indiceAles);
            fprintf('Construim mozaic ... %d biti de completat \n',p * q - l);
        end
    end
else
    nrTotalPiese = params.numarPieseMozaicOrizontala * params.numarPieseMozaicVerticala;
    nrPieseAdaugate = 0;
    medii = zeros(N);
    pieseGri = uint8(zeros(H, W, N));
    for i = 1:N
        pieseGri(:, :, i) = rgb2gray(params.pieseMozaic(:, :, :, i));
        medii(i) = mean2(pieseGri(:, :, i));
    end
    for l = 1:(p * q)
        i = pairs(perm(l), 1);
        j = pairs(perm(l), 2);
        if imgMozaic(i, j) == 0 || i == p || j == q
    
            medie = mean2(params.imgReferintaRedimensionata(i : i + H - 1,j : j + W - 1));
            distances = zeros(N);
            for k = 1:N
                distances(k) = abs(medie - medii(k));
            end
            [~ , indiceAles] = min(distances(:,1));
            imgMozaic(i : i + H - 1,j : j + W - 1) = pieseGri(:,:,indiceAles);
            fprintf('Construim mozaic ... %d biti de completat \n',p * q - l);
        end
    end
end
end

