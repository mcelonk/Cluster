format long g


XYZ = importdata('test_pointcloud.xyz');
T = 40; %velikost mezní hodnoty


D = pdist(XYZ,'euclidean'); %vzdálenost jednotlivých bodù
Dsquare = squareform(D); %vytvoøení matice
kmax = size(XYZ,1);

clusters = {}; %list clusterù
check_values = []%hodnoty vložené do clusteru

%vytvoøení clusteru podle mezní hodnoty T
for k = 1:kmax
    currentCluster = [];
    for i = 1:kmax
        pt = XYZ(k,:);
        pt_d = Dsquare(k,i); %vzdálenost konkrétních dvou bodù
        if pt_d < T %podmínka pro kontrolu
             if ~ismember(XYZ(i,:),check_values) %kontroluje jestli už není bod v clusteru
                currentCluster = [currentCluster; XYZ(i,:)];
                currentCluster(:,:,:);                 
             end
        end
    end
    if ~isempty(currentCluster) %pøidání clusteru do listu clusterù
        clusters = [ clusters ; {currentCluster} ];
        check_values = cell2mat(clusters);
    end
end

clusterCount = numel(clusters); %poèet clusterù
pointsNotInClusters = setdiff(XYZ,check_values,'rows'); %hodnoty, které nejsou v clusteru
isize = size(pointsNotInClusters,1);

%vložení zbylých bodù do clusteru podle vzdálenosti k centroidu
for i = 1:isize 
    dist_tab = [];
    for k = 1:clusterCount %každý cluster vytvoøený výše
           if size(clusters{k}(:,1),1) < 2 %kontrola pokud je cluster tvoøen pouze jedním bodem
               X = clusters{k,1}(:,:);
           else
               X = mean(clusters{k,1}(:,:)); %centroid
           end
           Y = pointsNotInClusters(i,:);
           XY = [X;Y];
           cdist = pdist(XY,'euclidean'); %vzdálenost k centroidu
           dist_tab = [dist_tab;cdist]; %tabulka vzdáleností
    end
    %pøiøazení do daného clusteru
    [M,I] = min(dist_tab);
    clusters{I,1} = [clusters{I,1};pointsNotInClusters(i,:)];
end



fprintf(1,'Number of clusters: %lu\n',clusterCount);

%Vykreslení clusterù
cc=hsv(clusterCount); % create different colour codes for every cluster
cc = cc(randperm(clusterCount),:); % randomise the colour codes so that neighbouring clusters don't have too similarly looking colours
h1 = figure('Name','Clusters');
hold on;
%scatter(XYZ(:,1),XYZ(:,2),20,'filled','o','CData',[.8,.8,.8]); originální
%body šeïì
for k=1:clusterCount
    plot(clusters{k,1}(:,1),clusters{k,1}(:,2),'o','Color',cc(k,:),'MarkerFaceColor',cc(k,:));
end
axis equal;

%writematrix(---); funkce pro uložení 