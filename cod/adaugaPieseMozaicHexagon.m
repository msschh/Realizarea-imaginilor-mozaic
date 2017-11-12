function imgMozaic = adaugaPieseMozaicHexagon( params )

[H,W,~,N] = size(params.pieseMozaic);
[h,w,c] = size(params.imgReferintaRedimensionata);
imgMozaic = uint8(zeros(h, w, c));

filter = fspecial('gaussian', 2 * W + 1, 1);
imagineExtinsa = imfilter(params.imgReferintaRedimensionata, filter, 'symmetric', 'full');
imagineExtinsa(H / 2 + 1 : h + H / 2, H / 2 + 1 : w + H / 2, :) = params.imgReferintaRedimensionata;
matriceHexagon = uint8(zeros(H, W));
for i=1:H
    for j=1:W
        if i > (H/2 - j) && i < (H/2 + j)
            matriceHexagon(i, j) = 1;
        end
    end
end
for i=1:H
    for j=1:W/2
        matriceHexagon(i, W - j + 1) = matriceHexagon(i, j);
    end
end

switch(params.criteriu)
    case 'distantaCuloareMedie'
        if c == 3
            medii = zeros(N, 3);
            for i = 1:N
                medii(i, 1) = mean2(matriceHexagon .* params.pieseMozaic(:, :, 1, i));
                medii(i, 2) = mean2(matriceHexagon .* params.pieseMozaic(:, :, 2, i));
                medii(i, 3) = mean2(matriceHexagon .* params.pieseMozaic(:, :, 3, i));
            end
            nr = 0;
            i = 0;
            while i < h + 2 * H
                if mod(nr,2) == 0
                    j = 0;
                else
                    j = W - floor(W/3) - 1;
                end
                while j < w + 3 * H / 2
                    mR = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W, 1));
                    mG = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W, 2));
                    mB = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W, 3));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = (mR - medii(k,1)).^2 + (mG - medii(k,2)).^2 + (mB - medii(k,3)).^2;
                    end
                    [~ , indiceAles] = min(distances(:,1));
                    imagineExtinsa(i+1:i+H, j+1:j+W, 1) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W, 1) + matriceHexagon .* params.pieseMozaic(:,:,1,indiceAles);
                    imagineExtinsa(i+1:i+H, j+1:j+W, 2) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W, 2) + matriceHexagon .* params.pieseMozaic(:,:,2,indiceAles); 
                    imagineExtinsa(i+1:i+H, j+1:j+W, 3) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W, 3) + matriceHexagon .* params.pieseMozaic(:,:,3,indiceAles);  
                    j = j + W + W - H;
                end
                i = i + H / 2;
                nr = nr + 1;
            end
        else
            medii = zeros(N);
            pieseGri = uint8(zeros(H, W, N));
            for i = 1:N
                pieseGri(:, :, i) = rgb2gray(params.pieseMozaic(:, :, :, i));
                medii(i) = mean2(matriceHexagon .* pieseGri(:, :, i));
            end
            nr = 0;
            i = 0;
            while i < h + 2 * H
                if mod(nr,2) == 0
                    j = 0;
                else
                    j = W - floor(W/3) - 1;
                end
                while j < w + 3 * H / 2
                    medie = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = abs(medie - medii(k));
                    end
                    [~ , indiceAles] = min(distances(:,1));
                    imagineExtinsa(i+1:i+H, j+1:j+W) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W) + matriceHexagon .* pieseGri(:,:,indiceAles);
                    j = j + W + W - H;
                end
                i = i + H / 2;
                nr = nr + 1;
            end
        end
        imgMozaic = imagineExtinsa(H + 1: h + H, H + 1: w + H, :);
    case 'distantaCuloareMedieDiferite'
        if c == 3
            medii = zeros(N, 3);
            matriceAparitii = zeros(floor(2.5 * params.numarPieseMozaicVerticala), floor(2.5 * params.numarPieseMozaicOrizontala));
            for i = 1:N
                medii(i, 1) = mean2(matriceHexagon .* params.pieseMozaic(:, :, 1, i));
                medii(i, 2) = mean2(matriceHexagon .* params.pieseMozaic(:, :, 2, i));
                medii(i, 3) = mean2(matriceHexagon .* params.pieseMozaic(:, :, 3, i));
            end
            nr = 0;
            i = 0;
            while i < h + 2 * H
                if mod(nr,2) == 0
                    j = 0;
                    l = 1;
                else
                    j = W - floor(W/3) - 1;
                    l = 2;
                end
                while j < w + 3 * H / 2
                    mR = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W, 1));
                    mG = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W, 2));
                    mB = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W, 3));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = (mR - medii(k,1)).^2 + (mG - medii(k,2)).^2 + (mB - medii(k,3)).^2;
                    end
                    [~, ind] = sort(distances,'ascend');
                    indiceAles = ind(1);
                    if nr == 0
                    elseif nr == 1
                        if matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(2);
                            if matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                                indiceAles = ind(3);
                            end
                        end
                    elseif j == 0
                        if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(2);
                            if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                                indiceAles = ind(3);
                            end
                        end
                    elseif matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(2);
                        if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(3);
                            if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                                indiceAles = ind(4);
                            end
                        end
                    end
                    matriceAparitii(nr + 1, l) = indiceAles;
                    l = l + 2;
                    imagineExtinsa(i+1:i+H, j+1:j+W, 1) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W, 1) + matriceHexagon .* params.pieseMozaic(:,:,1,indiceAles);
                    imagineExtinsa(i+1:i+H, j+1:j+W, 2) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W, 2) + matriceHexagon .* params.pieseMozaic(:,:,2,indiceAles); 
                    imagineExtinsa(i+1:i+H, j+1:j+W, 3) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W, 3) + matriceHexagon .* params.pieseMozaic(:,:,3,indiceAles);  
                    j = j + W + W - H;
                end
                i = i + H / 2;
                nr = nr + 1;
            end
        else
            medii = zeros(N);
            pieseGri = uint8(zeros(H, W, N));
            for i = 1:N
                pieseGri(:, :, i) = rgb2gray(params.pieseMozaic(:, :, :, i));
                medii(i) = mean2(matriceHexagon .* pieseGri(:, :, i));
            end
            matriceAparitii = zeros(floor(2.5 * params.numarPieseMozaicVerticala), floor(2.5 * params.numarPieseMozaicOrizontala));
            nr = 0;
            i = 0;
            while i < h + 2 * H
                if mod(nr,2) == 0
                    j = 0;
                    l = 1;
                else
                    j = W - floor(W/3) - 1;
                    l = 2;
                end
                while j < w + 3 * H / 2
                    medie = mean2(matriceHexagon .* imagineExtinsa(i + 1 : i + H,j + 1 : j + W));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = abs(medie - medii(k));
                    end
                    
                    [~, ind] = sort(distances,'ascend');
                    indiceAles = ind(1);
                    if nr == 0
                    elseif nr == 1
                        if matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(2);
                            if matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                                indiceAles = ind(3);
                            end
                        end
                    elseif j == 0
                        if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(2);
                            if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                                indiceAles = ind(3);
                            end
                        end
                    elseif matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(2);
                        if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                            indiceAles = ind(3);
                            if matriceAparitii(nr - 1, l) == indiceAles || matriceAparitii(nr, l - 1) == indiceAles || matriceAparitii(nr, l + 1) == indiceAles
                                indiceAles = ind(4);
                            end
                        end
                    end
                    matriceAparitii(nr + 1, l) = indiceAles;
                    l = l + 2;
                    imagineExtinsa(i+1:i+H, j+1:j+W) = (1 - matriceHexagon) .* imagineExtinsa(i+1:i+H, j+1:j+W) + matriceHexagon .* pieseGri(:,:,indiceAles);
                    j = j + W + W - H;
                end
                i = i + H / 2;
                nr = nr + 1;
            end
        end
        imgMozaic = imagineExtinsa(H + 1: h + H, H + 1: w + H, :);
end

