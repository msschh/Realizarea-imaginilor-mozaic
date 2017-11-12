function params = calculeazaDimensiuniMozaic(params)
%calculeaza dimensiunile mozaicului
%obtine si imaginea de referinta redimensionata avand aceleasi dimensiuni
%ca mozaicul

%completati codul Matlab
...
[h, w, ~] = size(params.pieseMozaic);
[H, W, ~] = size(params.imgReferinta);

%calculeaza automat numarul de piese pe verticala
params.numarPieseMozaicVerticala = floor(params.numarPieseMozaicOrizontala * w * H / h / W); 
params.imgReferintaRedimensionata = imresize(params.imgReferinta, [(params.numarPieseMozaicVerticala * h) (w * params.numarPieseMozaicOrizontala)]);