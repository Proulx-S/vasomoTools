function volTr = flashVolTr(rfTR,matSz,PAT,refLine)
% seperated GRAPPA ACS: refLine = 0
% integrated GRAPPA ACS: refLine = n lines
volTr = (refLine+(matSz-refLine)/PAT)*rfTR;
