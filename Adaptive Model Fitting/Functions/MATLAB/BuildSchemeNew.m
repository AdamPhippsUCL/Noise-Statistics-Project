
function scheme = BuildSchemeNew(scans)
% Build a scheme structure
% scans(indx, delta, Delta, bval)


for scanIndx = 1:size(scans,1)

    delta = scans(scanIndx, 1);
    Delta = scans(scanIndx,2);
    bval = scans(scanIndx,3);
    TE = scans(scanIndx, 4);
    TR = scans(scanIndx, 5);
    NSA = scans(scanIndx, 6);
    Rav = scans(scanIndx, 7);

    scheme(scanIndx).delta = delta;
    scheme(scanIndx).DELTA = Delta ;
    scheme(scanIndx).bval = bval;
    scheme(scanIndx).G = stejskal(delta,Delta,bval=bval);
    scheme(scanIndx).TE = TE;
    scheme(scanIndx).TR = TR;
    scheme(scanIndx).NSA = NSA;
    scheme(scanIndx).Rav = Rav;
end
    
end