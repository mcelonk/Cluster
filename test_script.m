format long g


XYZ = importdata('test_pointcloud.xyz');
T = 40; %velikost mezn� hodnoty


D = pdist(XYZ,'euclidean'); %vzd�lenost jednotliv�ch bod�
Dsquare = squareform(D); %vytvo�en� matice
kmax = size(XYZ,1);

clusters = {}; %list cluster�
check_values = []%hodnoty vlo�en� do clusteru

%vytvo�en� clusteru podle mezn� hodnoty T
for k = 1:kmax
    currentCluster = [];
    for i = 1:kmax
        pt = XYZ(k,:);
        pt_d = Dsquare(k,i); %vzd�lenost konkr�tn�ch dvou bod�
        if pt_d < T %podm�nka pro kontrolu
             if ~ismember(XYZ(i,:),check_values) %kontroluje jestli u� nen� bod v clusteru
                currentCluster = [currentCluster; XYZ(i,:)];
                currentCluster(:,:,:);                 
             end
        end
    end
    if ~isempty(currentCluster) %p�id�n� clusteru do listu cluster�
        clusters = [ clusters ; {currentCluster} ];
        check_values = cell2mat(clusters);
    end
end

clusterCount = numel(clusters); %po�et cluster�
pointsNotInClusters = setdiff(XYZ,check_values,'rows'); %hodnoty, kter� nejsou v clusteru
isize = size(pointsNotInClusters,1);

%vlo�en� zbyl�ch bod� do clusteru podle vzd�lenosti k centroidu
for i = 1:isize 
    dist_tab = [];
    for k = 1:clusterCount %ka�d� cluster vytvo�en� v��e
           if size(clusters{k}(:,1),1) < 2 %kontrola pokud je cluster tvo�en pouze jedn�m bodem
               X = clusters{k,1}(:,:);
           else
               X = mean(clusters{k,1}(:,:)); %centroid
           end
           Y = pointsNotInClusters(i,:);
           XY = [X;Y];
           cdist = pdist(XY,'euclidean'); %vzd�lenost k centroidu
           dist_tab = [dist_tab;cdist]; %tabulka vzd�lenost�
    end
    %p�i�azen� do dan�ho clusteru
    [M,I] = min(dist_tab);
    clusters{I,1} = [clusters{I,1};pointsNotInClusters(i,:)];
end



fprintf(1,'Number of clusters: %lu\n',clusterCount);

%Vykreslen� cluster�
cc=hsv(clusterCount); % create different colour codes for every cluster
cc = cc(randperm(clusterCount),:); % randomise the colour codes so that neighbouring clusters don't have too similarly looking colours
h1 = figure('Name','Clusters');
hold on;
%scatter(XYZ(:,1),XYZ(:,2),20,'filled','o','CData',[.8,.8,.8]); origin�ln�
%body �e��
for k=1:clusterCount
    plot(clusters{k,1}(:,1),clusters{k,1}(:,2),'o','Color',cc(k,:),'MarkerFaceColor',cc(k,:));
end
axis equal;

%writematrix(---); funkce pro ulo�en� 