function roiGrpMrgd = mergeRoi(roiGrp,commonFrameFlag)
    if ~exist('commonFrameFlag','var'); commonFrameFlag = []; end
    if isempty(commonFrameFlag);        commonFrameFlag = 1 ; end
    
    
    % Extract full-size underlay image from roi
    UL = [roiGrp{:}]; UL = [UL.im]; UL = [UL.base];
    UL = unique({UL.fName});
    if length(UL) ~= 1; dbstack; error('UL should be a single file name'); end;
    UL = char(UL);

    % use common frame
    if commonFrameFlag
        base = [roiGrp{:}]; base = [base.im]; base = [base.base];
        base(end+1).fName = unique({base.fName}); if length(base(end).fName) ~= 1; dbstack; error('multiple fNames in roiGrp'); end; base(end).fName = char(base(end).fName);
        base(end  ).x  = cat(1,base.x)    ; base(end).x  = [min(base(end).x(:,1)) max(base(end).x(:,2))];
        base(end  ).y  = cat(1,base.y)    ; base(end).y  = [min(base(end).y(:,1)) max(base(end).y(:,2))];
        maxSize = max(range(base(end).x),range(base(end).y)); base(end).x = round(mean(base(end).x) + [-maxSize/2 maxSize/2]); base(end).y = round(mean(base(end).y) + [-maxSize/2 maxSize/2]);
        base(end  ).im = MRIread(char(UL)); base(end).im = base(end).im.vol(base(end).y(1):base(end).y(2),base(end).x(1):base(end).x(2),:,:);
    end

    roiGrpMrgd = cell(size(roiGrp));
    for g = 1:length(roiGrp)
        if isempty(roiGrp{g}); continue; end

        % use group-specific frames
        if ~commonFrameFlag
            base = [roiGrp{g}]; base = [base.im]; base = [base.base];
            base(end+1).x  = cat(1,base.x)    ; base(end).x  = [min(base(end).x(:,1)) max(base(end).x(:,2))];
            base(end  ).y  = cat(1,base.y)    ; base(end).y  = [min(base(end).y(:,1)) max(base(end).y(:,2))];
            maxSize = max(range(base(end).x),range(base(end).y)); base(end).x = round(mean(base(end).x) + [-maxSize/2 maxSize/2]); base(end).y = round(mean(base(end).y) + [-maxSize/2 maxSize/2]);
            base(end  ).im = MRIread(char(UL)); base(end).im = base(end).im.vol(base(end).y(1):base(end).y(2),base(end).x(1):base(end).x(2),:,:);
        end


        % merge roi fields
        roiGrpMrgd{g}.poly  = [roiGrp{g}.poly];
        roiGrpMrgd{g}.mask  = any(cat(4,roiGrp{g}.mask),4);
        roiGrpMrgd{g}.class = unique({roiGrp{g}.class});
        if length(roiGrpMrgd{g}.class) == 1
            roiGrpMrgd{g}.class = char(roiGrpMrgd{g}.class);
        else
            dbstack; error('multiple classes in roiGrp');
        end
        roiGrpMrgd{g}.id = inf;
        roiGrpMrgd{g}.label = [roiGrpMrgd{g}.class 'All'];
        roiGrpMrgd{g}.com = {roiGrp{g}.com};
        roiGrpMrgd{g}.im.base = base(end);


        % merge data fields
        fieldList = fields(roiGrp{g}(1).vec.mt);
        for fi = 1:length(fieldList)
            sz = size(roiGrp{g}(1).vec.mt.(fieldList{fi}).vec,1:8); sz(6) = length(roiGrp{g});
            roiGrpMrgd{g}.vec.mt.(fieldList{fi}).vec = zeros(sz);
            for r = 1:length(roiGrp{g})
                roiGrpMrgd{g}.vec.mt.(fieldList{fi}).vec(:,:,:,:,:,r,:,:) = mean(roiGrp{g}(r).vec.mt.(fieldList{fi}).vec,6);
            end

            % copy info fields
            fieldList2 = fields(roiGrp{g}(1).vec.mt.(fieldList{fi}));
            fieldList2(ismember(fieldList2,{'vec' 'vecAv' 'vecEr'})) = [];
            for ffi = 1:length(fieldList2)
                roiGrpMrgd{g}.vec.mt.(fieldList{fi}).(fieldList2{ffi}) = roiGrp{g}(1).vec.mt.(fieldList{fi}).(fieldList2{ffi});
            end

            roiGrpMrgd{g}.vec.mt.psd.info = strsplit(roiGrpMrgd{g}.vec.mt.psd.info,' x ');
            roiGrpMrgd{g}.vec.mt.psd.info{6} = 'roi';
            roiGrpMrgd{g}.vec.mt.psd.info = strjoin(roiGrpMrgd{g}.vec.mt.psd.info,' x ');
        end
    end


    