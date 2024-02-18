function [archive_final] = select_long_pairs (archive)

ori_pref = [];
Is_resp = [];
dff_prefori_S0 = [];
dff_prefori_S1 = [];
longID = [];

ori_pref (1,:) = archive {3,1};
ori_pref (2,:) = archive {3,2};

Is_resp (1,:) = archive {5,1};
Is_resp (2,:) = archive {5,2};

dff_prefori_S0 (:,:) = archive {4,1};
dff_prefori_S1 (:,:) = archive {4,2};

longID (1,:) = archive {1,1};
longID (2,:) = archive {1,2};

%% Eliminate pairs that differ substantially on their preferred
% orientation

col = find(abs(ori_pref(1,:) - ori_pref(2,:))>1 & abs(ori_pref(1,:) - ori_pref(2,:))<5);

ori_pref (:,col) = [];
Is_resp (:,col) = [];
dff_prefori_S0 (:,col) = [];
dff_prefori_S1 (:,col) = [];
longID (:,col)=[];

clear col

%% Eliminate pairs that are not responsive in any of the sessions

col = find(Is_resp (1,:)+Is_resp (2,:) ==0);

ori_pref (:,col) = [];
Is_resp (:,col) = [];
dff_prefori_S0 (:,col) = [];
dff_prefori_S1 (:,col) = [];
longID (:,col)=[];

clear col
%% Eliminate pairs that are not responsive in S0 and S1

col = abs(Is_resp(1,:)-Is_resp(2,:)) ==(1);

ori_pref (:,col) = [];
Is_resp (:,col) = [];
dff_prefori_S0 (:,col) = [];
dff_prefori_S1 (:,col) = [];
longID (:,col)=[];

clear col

%% selected pairs
archive_final {1,1} = longID (1,:);
archive_final {1,2} = longID (2,:);
archive_final {2,1} = Is_resp (1,:);
archive_final {2,2} = Is_resp (2,:);
archive_final {3,1} = dff_prefori_S0;
archive_final {3,2} = dff_prefori_S1;
archive_final {4,1} = ori_pref (1,:);
archive_final {4,2} = ori_pref (2,:);




