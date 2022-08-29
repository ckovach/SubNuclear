

%%% Project view of contacts onto a plane in front of any obstructions.
cp = get(gca,'cameraposition');

dct = get(gca,'cameratarget')-cp;
uct = dct./sqrt(sum(dct.^2));

pmat = eye(3)- uct'*uct;

 %dcam = preop.transforms(3).tr(prx)-repmat(cp,size(prx,1),1);
dcam = prx(geti,:)-repmat(cp,size(prx(geti,:),1),1);
md = mean(sqrt(sum(dcam.^2,2)));
pflat = dcam*pmat + repmat(cp+.8*md*uct,size(dcam,1),1);




