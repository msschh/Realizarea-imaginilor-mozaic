function imgMozaic = adaugaPieseMozaicPeCaroiaj(params)
%
%tratati si cazul in care imaginea de referinta este gri (are numai un canal)

imgMozaic = uint8(zeros(size(params.imgReferintaRedimensionata)));
[H,W,~,N] = size(params.pieseMozaic);
[~,~,c] = size(params.imgReferintaRedimensionata);

switch(params.criteriu)
    case 'aleator'
        %pune o piese aleatoare in mozaic, nu tine cont de nimic
        nrTotalPiese = params.numarPieseMozaicOrizontala * params.numarPieseMozaicVerticala;
        nrPieseAdaugate = 0;
        for i =1:params.numarPieseMozaicVerticala
            for j=1:params.numarPieseMozaicOrizontala
                %alege un indice aleator din cele N
                indice = randi(N);    
                imgMozaic((i-1)*H+1:i*H,(j-1)*W+1:j*W,:) = params.pieseMozaic(:,:,:,indice);
                nrPieseAdaugate = nrPieseAdaugate+1;
                fprintf('Construim mozaic ... %2.2f%% \n',100*nrPieseAdaugate/nrTotalPiese);
            end
        end
        
    case 'distantaCuloareMedie'
        if c == 3
            nrTotalPiese = params.numarPieseMozaicOrizontala * params.numarPieseMozaicVerticala;
            nrPieseAdaugate = 0;
            medii = zeros(N, 3);
            for i = 1:N
                medii(i, 1) = mean2(params.pieseMozaic(:, :, 1, i));
                medii(i, 2) = mean2(params.pieseMozaic(:, :, 2, i));
                medii(i, 3) = mean2(params.pieseMozaic(:, :, 3, i));
            end
            for i = 1:params.numarPieseMozaicVerticala
                for j= 1:params.numarPieseMozaicOrizontala
                    mR = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W,1));
                    mG = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W,2));
                    mB = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W,3));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = (mR - medii(k,1)).^2 + (mG - medii(k,2)).^2 + (mB - medii(k,3)).^2;
                    end
                    [~ , indiceAles] = min(distances(:,1));
                    imgMozaic((i-1)*H+1:i*H,(j-1)*W+1:j*W,:) = params.pieseMozaic(:,:,:,indiceAles);
                    nrPieseAdaugate = nrPieseAdaugate+1;
                    fprintf('Construim mozaic ... %2.2f%% \n',100*nrPieseAdaugate/nrTotalPiese);
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
            for i = 1:params.numarPieseMozaicVerticala
                for j= 1:params.numarPieseMozaicOrizontala
                    medie = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = abs(medie - medii(k));
                    end
                    [~ , indiceAles] = min(distances(:,1));
                    imgMozaic((i-1)*H+1:i*H,(j-1)*W+1:j*W) = pieseGri(:,:,indiceAles);
                    nrPieseAdaugate = nrPieseAdaugate+1;
                    fprintf('Construim mozaic ... %2.2f%% \n',100*nrPieseAdaugate/nrTotalPiese);
                end
            end
        end
    case 'distantaCuloareMedieDiferite'
        if c == 3
            nrTotalPiese = params.numarPieseMozaicOrizontala * params.numarPieseMozaicVerticala;
            nrPieseAdaugate = 0;
            medii = zeros(N, 3);
            matriceAparitii = zeros(params.numarPieseMozaicVerticala, params.numarPieseMozaicOrizontala);
            for i = 1:N
                medii(i, 1) = mean2(params.pieseMozaic(:, :, 1, i));
                medii(i, 2) = mean2(params.pieseMozaic(:, :, 2, i));
                medii(i, 3) = mean2(params.pieseMozaic(:, :, 3, i));
            end
            for i = 1:params.numarPieseMozaicVerticala
                for j= 1:params.numarPieseMozaicOrizontala
                    mR = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W,1));
                    mG = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W,2));
                    mB = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W,3));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = (mR - medii(k,1)).^2 + (mG - medii(k,2)).^2 + (mB - medii(k,3)).^2;
                    end
                    [~, ind] = sort(distances,'ascend');
                    indiceAles = ind(1);
                    if i == 1
                        if j ~= 1 && matriceAparitii(i, j - 1) == indiceAles
                            indiceAles = ind(2);
                        end
                    elseif j == 1
                        if matriceAparitii(i - 1, j) == indiceAles
                            indiceAles = ind(2);
                        end
                    elseif matriceAparitii(i, j - 1) == indiceAles || matriceAparitii(i - 1, j) == indiceAles
                            indiceAles = ind(2);
                        if matriceAparitii(i, j - 1) == indiceAles || matriceAparitii(i - 1, j) == indiceAles
                            indiceAles = ind(3);
                        end
                    end
                    matriceAparitii(i, j) = indiceAles;
                    imgMozaic((i-1)*H+1:i*H,(j-1)*W+1:j*W,:) = params.pieseMozaic(:,:,:,indiceAles);
                    nrPieseAdaugate = nrPieseAdaugate+1;
                    fprintf('Construim mozaic ... %2.2f%% \n',100*nrPieseAdaugate/nrTotalPiese);
                end
            end
        else
            nrTotalPiese = params.numarPieseMozaicOrizontala * params.numarPieseMozaicVerticala;
            nrPieseAdaugate = 0;
            medii = zeros(N);
            pieseGri = uint8(zeros(H, W, N));
            matriceAparitii = zeros(params.numarPieseMozaicVerticala, params.numarPieseMozaicOrizontala);
            for i = 1:N
                pieseGri(:, :, i) = rgb2gray(params.pieseMozaic(:, :, :, i));
                medii(i) = mean2(pieseGri(:, :, i));
            end
            for i = 1:params.numarPieseMozaicVerticala
                for j= 1:params.numarPieseMozaicOrizontala
                    medie = mean2(params.imgReferintaRedimensionata((i-1)*H+1:i*H,(j-1)*W+1:j*W));
                    distances = zeros(N);
                    for k = 1:N
                        distances(k) = abs(medie - medii(k));
                    end
                    
                    [~, ind] = sort(distances,'ascend');
                    indiceAles = ind(1);
                    if i == 1
                        if j ~= 1 && matriceAparitii(i, j - 1) == indiceAles
                            indiceAles = ind(2);
                        end
                    elseif j == 1
                        if matriceAparitii(i - 1, j) == indiceAles
                            indiceAles = ind(2);
                        end
                    elseif matriceAparitii(i, j - 1) == indiceAles || matriceAparitii(i - 1, j) == indiceAles
                            indiceAles = ind(2);
                        if matriceAparitii(i, j - 1) == indiceAles || matriceAparitii(i - 1, j) == indiceAles
                            indiceAles = ind(3);
                        end
                    end
                    matriceAparitii(i, j) = indiceAles;
                    imgMozaic((i-1)*H+1:i*H,(j-1)*W+1:j*W) = pieseGri(:,:,indiceAles);
                    nrPieseAdaugate = nrPieseAdaugate+1;
                    fprintf('Construim mozaic ... %2.2f%% \n',100*nrPieseAdaugate/nrTotalPiese);
                end
            end
        end
        
    otherwise
        printf('EROARE, optiune necunoscuta \n');
    
end
    
    
    
    
    
