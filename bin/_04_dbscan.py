import pandas as pd
import numpy as np
import os
import re
from sklearn.neighbors import NearestNeighbors
import glob

#os.chdir("/data/dragon226/awatson/dej_rip/bin/")

files = sorted(glob.glob("../output/dist_mats/fullmatrix_*"))

for csv in files:
    mega_mat = pd.read_csv(csv)
    mega_mat = mega_mat.to_numpy()
    print(mega_mat.shape)

    mds_nm = re.sub("dist_mats/fullmatrix", "mds_dists/vis_mds", csv)
    mds_dist = pd.read_csv(mds_nm)
    print(mds_dist.shape)
    def plot_k(k):
        knn_dist = NearestNeighbors(n_neighbors=k).fit(mega_mat.T)
        distances, _ = knn_dist.kneighbors(mega_mat.T, n_neighbors=k)
        distances = distances.reshape(-1)
        knn_dist_sorted = np.sort(distances)[::-1]
        n = len(knn_dist_sorted)
        points_mat = np.column_stack((np.arange(1, n + 1), knn_dist_sorted))
        line_vec = points_mat[-1, :] - points_mat[0, :]

        def point_to_line_dist(point, line_vec, line_point=points_mat[0, :]):
            return np.abs((line_vec[1] * (line_point[0] - point[0])) - ((line_point[1] - point[1]) * line_vec[0])) / np.sqrt(np.sum(line_vec**2))

        distances = np.apply_along_axis(point_to_line_dist, 1, points_mat, line_vec=line_vec)
        elbow_index = np.argmax(distances)
        elbow_point = round(knn_dist_sorted[elbow_index], 5)
        print("The 'elbow' point is at a distance of:", elbow_point)
        return elbow_point

    elbow_pts = [plot_k(k) for k in range(2, 6)]
    from sklearn.cluster import DBSCAN

    def sklearn_dbscan(X, eps, min_samples):
        # DBSCAN clustering
        db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=20).fit(X)
        labels = db.labels_
        return labels

    labels = []

    for num, eps in enumerate(elbow_pts):
        labels.append(sklearn_dbscan(mega_mat, eps, round(
            ########
            np.sqrt(len(mega_mat)#/2 # Tweak this to get a balance
                    )
                                                          )))

    df = pd.DataFrame(labels).T
    df = df + 1
    df.columns = ["K2", "K3", "K4", "K5"]
    df = pd.concat([mds_dist, df], axis=1)
    nm = csv.replace("dist_mats/fullmatrix_", "knn_dbs/eps_knn_")
    dir_name = os.path.dirname(nm)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    df.to_csv(nm, index=False)

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))  # Create a 2x2 grid of subplots

    columns = ["K2", "K3", "K4", "K5"]

    for ax, column in zip(axs.flatten(), columns):
        labels = df[column].unique()
        colors = plt.cm.viridis(np.linspace(0, 1, len(labels)))
        for label, color in zip(labels, colors):
            index = df[column] == label
            ax.scatter(df.loc[index, 'Coord1'], df.loc[index, 'Coord2'], color=color, label=f'{label}: {sum(index)}',  s=2)
        ax.set_title(column)
        ax.set_xlabel('Coord1')
        ax.set_ylabel('Coord2')
        ax.legend(title='Label: Count')

    plt.tight_layout()  # Adjust subplot parameters to give specified padding
        
    # Save the plot
    plot_filename = os.path.splitext(nm)[0] + '.png'  # replace the .csv with .png in the filename
    plt.savefig(plot_filename)