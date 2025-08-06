function cluster_ids = fixedNumConnection(pos, num_connect, outlier_fac)
if nargin < 3
    outlier_fac = 5;
end
if num_connect == 1
    cluster_ids = (1:size(pos, 1))';
    return;
end
[X, Y] = meshgrid(1:length(pos), 1:length(pos));
distance = sqrt((pos(:, 1)'-pos(:, 1)).^2 + (pos(:, 2)'-pos(:, 2)).^2);
distance(eye(length(distance)) == 1) = nan;
X(eye(length(X)) == 1) = nan;
Y(eye(length(Y)) == 1) = nan;


nearest_dist = min(distance, [], 2);
mnearest_dist = median(nearest_dist);
max_iter = ceil(size(pos, 1)/num_connect);
cluster_ids = nan(max_iter, num_connect);
search_list = 1:size(pos, 1);
corner_dist = sqrt(sum(pos.^2, 2));
[~, corner_dist_ids] = sort(corner_dist);
search_id = corner_dist_ids(1);
ci = 1;
while length(search_list) >= num_connect && ci < max_iter
    [~, sids] = sort(distance(search_id, :), 'ascend');
    sids = sids(ismember(sids, search_list));
    sids = sids(1:(num_connect-1));
    if any(distance(search_id, sids) > mnearest_dist*outlier_fac)
        search_list(search_list == search_id) = [];
        corner_dist_ids(corner_dist_ids == search_id) = [];
    else
        cids = [search_id sids];
        if ismember(cids, cluster_ids(:))
            break
        end
        cluster_ids(ci, :) = cids;
        search_list(ismember(search_list, cids)) = [];
        corner_dist_ids(ismember(corner_dist_ids, cids)) = [];
    end
    ci = ci +1 ;
    search_id = corner_dist_ids(1);
end

cluster_ids(any(isnan(cluster_ids), 2), :) = [];