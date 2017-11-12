function params = incarcaPieseMozaic(params)
%citeste toate cele N piese folosite la mozaic din directorul corespunzator
%toate cele N imagini au aceeasi dimensiune H x W x C, unde:
%H = inaltime, W = latime, C = nr canale (C=1  gri, C=3 color)
%functia intoarce pieseMozaic = matrice H x W x C x N in params
%pieseMoziac(:,:,:,i) reprezinta piese numarul i 

fprintf('Incarcam piesele pentru mozaic din director \n');
%completati codul Matlab
filelist = dir(strcat(params.numeDirector,'*.',params.tipImagine));
img = imread([params.numeDirector filelist(1).name]);
[h, w, c] = size(img);
pieseMozaic = uint8(zeros(h, w, c, length(filelist)));
for idxImg = 1:length(filelist)
    imgName = [params.numeDirector filelist(idxImg).name];
    pieseMozaic(:, :, :, idxImg) = imread(imgName);
end
%pt cifrar
%{
nr = 0;
for idxImg = 1:length(filelist)
    if isempty(strfind(filelist(idxImg).name, params.categorieCifrar10)) == 0
        nr = nr + 1;
    end
end
pieseMozaic = uint8(zeros(h, w, c, nr));
nr = 0;
for idxImg = 1:length(filelist)
    if isempty(strfind(filelist(idxImg).name, params.categorieCifrar10)) == 0
        imgName = [params.numeDirector filelist(idxImg).name];
        nr = nr + 1;
        pieseMozaic(:, :, :, nr) = imread(imgName);
    end
end
%}

if params.afiseazaPieseMozaic
    %afiseaza primele 100 de piese ale mozaicului
    figure,
    title('Primele 100 de piese ale mozaicului sunt:');
    idxImg = 0;
    for i = 1:10
        for j = 1:10
            idxImg = idxImg + 1;
            subplot(10,10,idxImg);
            imshow(pieseMozaic(:,:,:,idxImg));
        end
    end
    drawnow;
    pause(2);
end

params.pieseMozaic = pieseMozaic;