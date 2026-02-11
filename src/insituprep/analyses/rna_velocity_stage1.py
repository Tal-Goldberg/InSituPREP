from __future__ import annotations

import os
import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# =========================
# Helpers (unchanged logic)
# =========================

def save_to_csv(file_name, header, row_names, values):
    df = pd.DataFrame(values, columns=header, index=row_names)
    df.to_csv(file_name, index=True)


def plot_precision_recall(csv_file: Path, out_path: Path):
    # Read data from CSV file
    expression_levels = np.arange(0, 0.061, 0.001)

    with open(csv_file, mode="r") as f:
        lines = f.readlines()
        # header = lines[0].strip().split(",")
        line_p = lines[1]
        line_r = lines[2]

        values_p = line_p.strip().split(",")[1:]
        values_r = line_r.strip().split(",")[1:]

        precision_values = [float(x) for x in values_p]
        recall_values = [float(x) for x in values_r]

    # Calculate F-score
    f_scores = [(2 * p * r) / (p + r) for p, r in zip(precision_values, recall_values)]

    max_f_score = max(f_scores)
    max_f_score_index = f_scores.index(max_f_score)
    max_f_score_expression_level = expression_levels[max_f_score_index]

    # Plot precision, recall, and F-score
    plt.figure(figsize=(8, 6), constrained_layout=True)
    plt.plot(expression_levels, precision_values, label="Precision", linewidth=5)
    plt.plot(expression_levels, recall_values, label="Recall", linewidth=5)
    plt.plot(expression_levels, f_scores, label="F-score", linewidth=5)
    plt.scatter(max_f_score_expression_level, max_f_score, color="red", s=120, zorder=10)

    plt.xlabel("Normalized expression levels", fontsize=20, labelpad=16)
    plt.ylabel("Precision / Recall / F-score", fontsize=20, labelpad=16)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=16)
    plt.grid(True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, bbox_inches="tight", dpi=300, pad_inches=0.2)
    plt.close()


def create_boolean_vector(filter_indices, vector):
    for index in filter_indices:
        vector[index] = 0
    return vector


def load_ground_truth_genes(txt_path: Path) -> List[str]:
    genes = []
    with open(txt_path, "r") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)
    return genes


# =========================
# Core classes (same logic)
# =========================

class CellsDict(dict):
    """
    Dict of cell name -> Cell
    """

    def __init__(
        self,
        df_pc: pd.DataFrame,
        spliced_counts_path: Path,
        unspliced_counts_path: Path,
        neighbors_csv_path: Path,
        first_run: bool,
        percentage: float,
    ):
        super().__init__()
        self.df = df_pc
        self.number_of_cells = self.df.shape[0]

        spliced_df = pd.read_csv(spliced_counts_path, index_col=0, dtype={0: str})
        unspliced_df = pd.read_csv(unspliced_counts_path, index_col=0, dtype={0: str})

        self.neighbors_csv_path = neighbors_csv_path
        self.first_run = first_run
        self.PERCENTAGE = percentage

        if self.first_run:
            self.neighbors_df = None
        else:
            self.neighbors_df = pd.read_csv(self.neighbors_csv_path)

        self.genes = []
        self.genes_name = []

        self.S = self.initiate_mRNA_expression(spliced_df, is_first=True)
        self.U = self.initiate_mRNA_expression(unspliced_df)

        self.expression_levels = np.arange(0, 0.061, 0.001)
        self.precision_recall_df = pd.DataFrame(columns=self.expression_levels.tolist())
        self.precision_recall_df = self.precision_recall_df.reindex(["precision", "recall"])

    def initiate_mRNA_expression(self, df: pd.DataFrame, is_first=False):
        cells_name = df.index.to_list()
        genes_name = df.columns.to_list()[1:]  # SAME as original
        values = df.values

        expression_dict = {}
        for i in range(self.number_of_cells):
            expression_dict.update(
                {cells_name[i]: {gene: value for gene, value in zip(genes_name, values[i, 1:])}}
            )

        if is_first:
            self.init_genes(genes_name)

        return expression_dict

    def init_genes(self, genes_name):
        self.genes_name = genes_name
        for name in genes_name:
            self.genes.append(Gene(name))

    def set_dict(self):
        cells_name = self.df.index.tolist()
        self.df["cell_type"] = self.df["cell_type"].astype(str)
        cells_type = self.df["cell_type"].tolist()

        self.cell_types = list(set(cells_type))

        pc1 = self.df["PC1"]
        pc2 = self.df["PC2"]
        pc3 = self.df["PC3"]

        for i in range(self.number_of_cells):
            cell_name = cells_name[i]
            self.update(
                {
                    cell_name: Cell(
                        name=cell_name,
                        S0=self.S[cell_name],
                        U0=self.U[cell_name],
                        cell_type=cells_type[i],
                        pca=np.array([pc1[i], pc2[i], pc3[i]]),
                    )
                }
            )

    def save_distances_to_csv(self):
        cells_neighbors = {}
        for current_cell in self.values():
            distances = []
            for cell in self.values():
                distance = np.linalg.norm(current_cell.pca - cell.pca, axis=0)
                distances.append((distance, cell))
            sorted_list = sorted(distances, key=lambda x: x[0])
            sorted_cells_names = [x[1].name for x in sorted_list]
            cells_neighbors.update({current_cell.name: sorted_cells_names})

        df = pd.DataFrame(cells_neighbors)
        self.neighbors_csv_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(self.neighbors_csv_path, index=False)
        self.neighbors_df = pd.read_csv(self.neighbors_csv_path)

    def calculate_gammas(self, l: float, K: int):
        gammas = []
        filter_list_false = []
        filter_list_final = []

        for gene in self.genes:
            points = gene.cells_coordinates.values()
            x_values = [point[0] for point in points]
            y_values = [point[1] for point in points]

            gene.filter = False
            filter_list_false.append(gene.filter)

            if max(x_values) < l and max(y_values) < l:
                gene.filter = True

            filter_list_final.append(gene.filter)

            # Find extreme points
            max_x, max_y = max(x_values), max(y_values)
            min_x, min_y = 0, 0

            num_points = len(list(points))
            num_points_percent = int(num_points * self.PERCENTAGE)

            # top
            sorted_points = sorted(
                list(gene.cells_coordinates.values()),
                key=lambda point: ((point[0] - max_x) ** 2 + (point[1] - max_y) ** 2) ** 0.5,
            )
            top_points = sorted_points[:num_points_percent]

            # bottom
            sorted_points = sorted(
                list(gene.cells_coordinates.values()),
                key=lambda point: ((point[0] - min_x) ** 2 + (point[1] - min_y) ** 2) ** 0.5,
            )
            bottom_points = sorted_points[:num_points_percent]

            points2 = top_points + bottom_points
            x2 = np.array([p[0] for p in points2]).reshape(-1, 1)
            y2 = np.array([p[1] for p in points2])

            reg = LinearRegression(fit_intercept=False)
            reg.fit(x2, y2)
            slope = reg.coef_

            gene.set_gamma(slope[0])
            gammas.append(slope[0])

        return gammas, filter_list_false, filter_list_final

    def set_gammas_for_cells(self, gammas):
        for cell in self.values():
            cell.set_gammas(gammas)


class Gene:
    def __init__(self, name):
        self.name = name
        self.gamma = None
        self.cells_coordinates = {}
        self.filter = False

    def set_gamma(self, slope):
        self.gamma = slope

    def get_gamma(self):
        return self.gamma

    def set_coordinates(self, cell_name, coordinates):
        self.cells_coordinates.update({cell_name: coordinates})


class Cell:
    def __init__(self, name, S0, U0, cell_type, pca=None):
        self.name = name
        self.S0 = S0
        self.total_S_counts = self.get_total_counts(self.S0)
        self.s0_norm = self.normalization(self.S0, self.total_S_counts)
        self.s_knn = {}

        self.U0 = U0
        self.total_U_counts = self.get_total_counts(self.U0)
        self.u0_norm = self.normalization(self.U0, self.total_U_counts)
        self.u_knn = {}

        self.cell_type = cell_type
        self.neighbors = None
        self.pca = pca
        self.gammas = None

    def set_gammas(self, gammas):
        self.gammas = gammas

    def get_st(self, cells: CellsDict, t: float, filtered_genes: List[str]):
        normal_v, filtered_v = self.get_v_calculation(cells, filtered_genes)

        norm_st = np.array(list(self.s0_norm.values())) + normal_v * t
        filtered_norm_st = np.array(list(self.s0_norm.values())) + filtered_v * t

        st = np.array(list(self.S0.values())) + normal_v * t
        filtered_st = np.array(list(self.S0.values())) + filtered_v * t

        return norm_st, filtered_norm_st, st, filtered_st, normal_v, filtered_v

    def get_v_calculation(self, cells: CellsDict, filtered_genes: List[str]):
        normal_v = np.array(list(self.u0_norm.values())) - self.gammas * np.array(list(self.s0_norm.values()))
        filter_indices = [index for index, gene in enumerate(cells.genes) if gene.name in filtered_genes]
        filtered_v = create_boolean_vector(filter_indices, normal_v.copy())
        return normal_v, filtered_v

    def initialize_neighbors(self, cells: CellsDict, K: int):
        self.neighbors = self.get_PCA_neighbors(cells, K)

    def get_total_counts(self, expression_dict):
        total = 0
        for x in expression_dict.values():
            total += x
        return total

    def normalization(self, expression_dict, total_counts):
        return {gene: expression / total_counts for gene, expression in expression_dict.items()}

    def substitute_S_by_knn(self):
        sum_vec = 0
        total_neighbor_S_counts = 0
        for neighbor in self.neighbors:
            sum_vec += np.array(list(neighbor.S0.values()))
            total_neighbor_S_counts += neighbor.total_S_counts

        sum_vec = sum_vec / total_neighbor_S_counts
        self.s_knn = {gene: sum_vec[i] for i, gene in enumerate(self.S0.keys())}

    def substitute_U_by_knn(self):
        sum_vec = 0
        total_neighbor_U_counts = 0
        for neighbor in self.neighbors:
            sum_vec += np.array(list(neighbor.U0.values()))
            total_neighbor_U_counts += neighbor.total_U_counts

        sum_vec = sum_vec / total_neighbor_U_counts
        self.u_knn = {gene: sum_vec[i] for i, gene in enumerate(self.U0.keys())}

    def get_PCA_neighbors(self, cells: CellsDict, K: int):
        names = list(cells.neighbors_df[str(self.name)].head(K))
        neighbors = []
        for i in range(K):
            neighbors.append(cells[names[i]])
        return neighbors


# =========================
# Precision/Recall (same logic)
# =========================

def precision_and_recall(cells: CellsDict, l: float, ground_truth: List[str]):
    tp = 0
    fp = 0
    predicted_genes = [gene.name for gene in cells.genes if gene.filter == False]
    for gene in predicted_genes:
        if gene in ground_truth:
            tp += 1
        else:
            fp += 1
    fn = len(ground_truth) - tp

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    cells.precision_recall_df[l] = [precision, recall]


# =========================
# Stage functions (CLI calls)
# =========================

def _build_cells(
    pc_origin_path: Path,
    spliced_counts_path: Path,
    unspliced_counts_path: Path,
    neighbors_csv_path: Path,
    first_run: bool,
    percentage: float,
) -> CellsDict:
    cells = CellsDict(
        df_pc=pd.read_csv(pc_origin_path, index_col=0, dtype={0: str}),
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        neighbors_csv_path=neighbors_csv_path,
        first_run=first_run,
        percentage=percentage,
    )
    cells.set_dict()

    print(f"[DATA] Loaded {cells.number_of_cells} cells")
    print(f"[DATA] Loaded {len(cells.genes)} genes")

    if first_run:
        print("[KNN] neighbors.csv not found -> computing PCA neighbors")
        cells.save_distances_to_csv()
        cells.first_run = False
    else:
        print("[KNN] Using existing neighbors.csv")

    return cells


def _compute_knn_coordinates(cells: CellsDict, K: int):
    for cell in cells.values():
        cell.initialize_neighbors(cells, K)
        cell.substitute_S_by_knn()
        cell.substitute_U_by_knn()

        for j, (x, y) in enumerate(zip(cell.s_knn.values(), cell.u_knn.values())):
            cells.genes[j].set_coordinates(cell, (x, y))


def stage1_gamma_grid(
    pc_origin_path: Path,
    spliced_counts_path: Path,
    unspliced_counts_path: Path,
    out_dir: Path,
    k_values: List[int],
    quantiles: List[float],
    l_fixed: float,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[stage1-a] Gamma grid computation")
    print(f"[stage1-a] Saving outputs to: {out_dir}")

    neighbors_csv_path = out_dir / "neighbors.csv"

    gamma_results: Dict[int, Dict[float, np.ndarray]] = {}

    first_run = not neighbors_csv_path.exists()

    for k in k_values:
        gamma_results[k] = {}
        for q in quantiles:
            print(f"[stage1-a] Running k={k}, quantile={q}")
            cells = _build_cells(
                pc_origin_path=pc_origin_path,
                spliced_counts_path=spliced_counts_path,
                unspliced_counts_path=unspliced_counts_path,
                neighbors_csv_path=neighbors_csv_path,
                first_run=first_run,
                percentage=q,
            )
            first_run = False

            _compute_knn_coordinates(cells, K=k)

            gammas, _, _ = cells.calculate_gammas(l_fixed, K=k)
            print(f"[stage1-a] Done k={k}, quantile={q} (gammas computed)")
            gamma_results[k][q] = np.array(gammas)

    # Build per-gene dict of dataframes (same structure)
    gamma_df_gene_dict: Dict[str, pd.DataFrame] = {}
    for i in range(len(list(cells.genes))):
        gene_name = list(cells.genes)[i].name
        gamma_df_gene = pd.DataFrame(index=quantiles, columns=k_values)
        for k in k_values:
            for q in quantiles:
                gamma_df_gene.at[q, k] = gamma_results[k][q][i]
        gamma_df_gene_dict[gene_name] = gamma_df_gene

    with open(out_dir / "gamma_per_gene_per_quantile_per_k_dict.pkl", "wb") as f:
        pickle.dump(gamma_df_gene_dict, f)

    print(f"[stage1-a] Saved: {out_dir / 'gamma_per_gene_per_quantile_per_k_dict.pkl'}")


def stage1_gamma_plots_all_genes(
    pc_origin_path: Path,
    spliced_counts_path: Path,
    unspliced_counts_path: Path,
    out_dir: Path,
    k: int,
    quantile: float,
    l_fixed: float,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    neighbors_csv_path = out_dir / "neighbors.csv"
    plots_dir = out_dir / "gamma_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    print("[stage1-b] Gamma portrait plots for ALL genes")
    print(f"[stage1-b] k={k}, quantile={quantile}")
    print(f"[stage1-b] Plots dir: {plots_dir}")


    first_run = not neighbors_csv_path.exists()

    cells = _build_cells(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        neighbors_csv_path=neighbors_csv_path,
        first_run=first_run,
        percentage=quantile,
    )
    _compute_knn_coordinates(cells, K=k)

    gammas, _, _ = cells.calculate_gammas(l_fixed, K=k)
    gammas = np.array(gammas)

    print(f"[stage1-b] Computed gammas for {len(gammas)} genes")

    # dynamic colors per dataset (general)
    cell_types = [c.cell_type for c in cells.values()]
    unique_types = sorted(list(set(cell_types)))
    cmap = plt.get_cmap("tab20" if len(unique_types) <= 20 else "hsv")
    colors_map = {ct: cmap(i / max(1, len(unique_types) - 1)) for i, ct in enumerate(unique_types)}
    c = [colors_map[cells[cell_name].cell_type] for cell_name in cells.keys()]

    for gene in cells.genes:
        points = gene.cells_coordinates.values()
        x_values = [p[0] for p in points]
        y_values = [p[1] for p in points]

        plt.figure(figsize=(8, 6))
        plt.scatter(x_values, y_values, c=c)

        legend_handles = []
        for ct in unique_types:
            legend_handles.append(
                plt.Line2D(
                    [0], [0], marker="o", color="w", label=str(ct),
                    markerfacecolor=colors_map[ct], markersize=10
                )
            )
        plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
        plt.subplots_adjust(right=0.7)

        # SAME gamma fit logic as original
        max_x, max_y = max(x_values), max(y_values)
        min_x, min_y = 0, 0
        num_points = len(list(points))
        num_points_percent = int(num_points * quantile)

        sorted_top = sorted(list(gene.cells_coordinates.values()),
                            key=lambda p: ((p[0] - max_x) ** 2 + (p[1] - max_y) ** 2) ** 0.5)
        sorted_bottom = sorted(list(gene.cells_coordinates.values()),
                               key=lambda p: ((p[0] - min_x) ** 2 + (p[1] - min_y) ** 2) ** 0.5)

        top_points = sorted_top[:num_points_percent]
        bottom_points = sorted_bottom[:num_points_percent]
        points2 = top_points + bottom_points

        x2 = np.array([p[0] for p in points2]).reshape(-1, 1)
        y2 = np.array([p[1] for p in points2])

        reg = LinearRegression(fit_intercept=False)
        reg.fit(x2, y2)
        slope = reg.coef_

        x_regression = np.linspace(0, np.max(x2), len(points2))
        y_regression = slope * x_regression
        plt.plot(x_regression, y_regression, color="black", linewidth=3)

        plt.xlabel("Normalized cytoplasmic (spliced) RNA", fontsize=14, labelpad=10)
        plt.ylabel("Normalized nuclear (unspliced) RNA", fontsize=14, labelpad=10)
        plt.title(f"{gene.name}", fontsize=16, pad=5)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))

        plt.savefig(plots_dir / f"{gene.name}_K{k}.tif", dpi=300, bbox_inches="tight")
        plt.close()

    print("[stage1-b] Finished gamma plots")


def stage1_find_best_l_fscore(
    pc_origin_path: Path,
    spliced_counts_path: Path,
    unspliced_counts_path: Path,
    out_dir: Path,
    k: int,
    quantile: float,
    ground_truth_txt: Path,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    neighbors_csv_path = out_dir / "neighbors.csv"
    first_run = not neighbors_csv_path.exists()

    ground_truth = load_ground_truth_genes(ground_truth_txt)

    print("[stage1-c] best_l via precision/recall F-score scan")
    print(f"[stage1-c] k={k}, quantile={quantile}")
    print(f"[stage1-c] ground_truth genes: {len(ground_truth)}")

    cells = _build_cells(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        neighbors_csv_path=neighbors_csv_path,
        first_run=first_run,
        percentage=quantile,
    )
    _compute_knn_coordinates(cells, K=k)

    print(f"[stage1-c] Scanning {len(cells.expression_levels)} expression thresholds")

    # loop on cells.expression_levels EXACTLY like original
    for l in cells.expression_levels:
        gammas, _, _ = cells.calculate_gammas(float(l), K=k)
        _ = np.array(gammas)
        precision_and_recall(cells, float(l), ground_truth)

    pr_csv = out_dir / "precision_recall_unspliced.csv"
    cells.precision_recall_df.to_csv(pr_csv, index=True)

    # plot + compute best_l from file (same as original plot_precision_recall behavior)
    plot_path = out_dir / "Precision_Recall_F-score_plot.tif"
    plot_precision_recall(pr_csv, plot_path)

    print(f"[stage1-c] Saved: {pr_csv}")
    print(f"[stage1-c] Saved: {plot_path}")

    # Return best_l as computed from the CSV (same logic as plot_precision_recall)
    expression_levels = np.arange(0, 0.061, 0.001)
    with open(pr_csv, "r") as f:
        lines = f.readlines()
        values_p = lines[1].strip().split(",")[1:]
        values_r = lines[2].strip().split(",")[1:]
        precision_values = [float(x) for x in values_p]
        recall_values = [float(x) for x in values_r]
    f_scores = [(2 * p * r) / (p + r) if (p + r) > 0 else 0.0 for p, r in zip(precision_values, recall_values)]
    best_idx = int(np.argmax(f_scores))
    best_l = float(expression_levels[best_idx])

    params_path = out_dir / "stage1_selected_params.txt"
    with open(params_path, "w") as f:
        f.write("# RNA velocity – selected parameters\n")
        f.write("# Generated by stage1-c\n\n")
        f.write(f"k_neighbors: {k}\n")
        f.write(f"quantile_outliers: {quantile}\n")
        f.write(f"l_normalized_expression_cutoff: {best_l}\n")

    print(f"[stage1-c] best_l = {best_l}")
    print(f"[stage1-c] Saved: {params_path}")

    return best_l


def stage1_t_preview(
    pc_origin_path: Path,
    spliced_counts_path: Path,
    unspliced_counts_path: Path,
    out_dir: Path,
    k: int,
    quantile: float,
    ground_truth_txt: Path,
    best_l: float,
    T: Optional[List[float]] = None,
    T2: Optional[List[float]] = None,
    n_cells: int = 500,
    cell_type: str = "randomly selected 500 cells",
):
    out_dir.mkdir(parents=True, exist_ok=True)
    neighbors_csv_path = out_dir / "neighbors.csv"
    first_run = not neighbors_csv_path.exists()

    ground_truth = load_ground_truth_genes(ground_truth_txt)

    def _tkey(t: float) -> float:
        return round(float(t), 1)  # 0.2 grid -> 1 decimal is enough

    cells = _build_cells(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        neighbors_csv_path=neighbors_csv_path,
        first_run=first_run,
        percentage=quantile,
    )
    _compute_knn_coordinates(cells, K=k)

    gammas, _, _ = cells.calculate_gammas(best_l, K=k)
    gammas = np.array(gammas)

    # final gene sets (same as original)
    genes_pass = [g.name for g in cells.genes if g.filter == False]
    not_filtered_genes = genes_pass
    filtered_genes = [g.name for g in cells.genes if g.name not in not_filtered_genes]

    # default time grids (same values)
    if T is None:
        T = [round(x, 1) for x in np.arange(0.2, 12, 0.2)]
    if T2 is None:
        T2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    T2 = [_tkey(x) for x in T2]

    print("[stage1-d] Computing s(t) across multiple t values")
    print(f"[stage1-d] T grid size: {len(T)}")
    print(f"[stage1-d] Arrow plots for t values: {T2}")


    change_in_sum_s_per_cell_all_t = {}
    all_gene_level_dfs_per_t = {}
    df_sum_s_dict = {}

    print("[stage1-d] Building gene-level DataFrames per t ...")

    for t in T:
        t_val = float(t)          # value used for numerical computation
        tk = _tkey(t_val)      # dictionary key (rounded to avoid float precision issues)

        # Print only every ~1.0 step (T is in 0.2 increments)
        if abs((tk * 10) % 10) < 1e-9:
            print(f"[stage1-d] t = {tk}")

        sum_s_dict = {}
        S0_dict = {}
        S0_norm_dict = {}
        st_dict = {}
        st_filtered_dict = {}
        st_norm_dict = {}
        st_norm_filtered_dict = {}
        normal_v_dict = {}
        filtered_v_dict = {}

        for cell in cells.values():
            cell.set_gammas(gammas)

            S0 = np.array(list(cell.S0.values()))
            S0_norm = np.array(list(cell.s0_norm.values()))

            # Use the raw float value for the actual s(t) computation
            norm_st, filtered_norm_st, st, filtered_st, normal_v, filtered_v = cell.get_st(
                cells, t_val, filtered_genes
            )

            sum_S0 = np.sum(S0)
            sum_S0_norm = np.sum(S0_norm)
            sum_st = np.sum(st)
            sum_st_filtered = np.sum(filtered_st)
            sum_st_norm = np.sum(norm_st)
            sum_st_norm_filtered = np.sum(filtered_norm_st)

            sum_s_dict[cell.name] = [
                sum_S0,
                sum_S0_norm,
                sum_st,
                sum_st_filtered,
                sum_st_norm,
                sum_st_norm_filtered,
            ]

            S0_dict[cell.name] = S0
            S0_norm_dict[cell.name] = S0_norm
            st_dict[cell.name] = st
            st_filtered_dict[cell.name] = filtered_st
            st_norm_dict[cell.name] = norm_st
            st_norm_filtered_dict[cell.name] = filtered_norm_st
            normal_v_dict[cell.name] = normal_v
            filtered_v_dict[cell.name] = filtered_v

        df_S0 = pd.DataFrame.from_dict(S0_dict, orient="index", columns=cells.genes_name)
        df_S0_norm = pd.DataFrame.from_dict(S0_norm_dict, orient="index", columns=cells.genes_name)
        df_st = pd.DataFrame.from_dict(st_dict, orient="index", columns=cells.genes_name)
        df_st_filtered = pd.DataFrame.from_dict(st_filtered_dict, orient="index", columns=cells.genes_name)
        df_st_norm = pd.DataFrame.from_dict(st_norm_dict, orient="index", columns=cells.genes_name)
        df_st_norm_filtered = pd.DataFrame.from_dict(st_norm_filtered_dict, orient="index", columns=cells.genes_name)
        df_normal_v = pd.DataFrame.from_dict(normal_v_dict, orient="index", columns=cells.genes_name)
        df_filtered_v = pd.DataFrame.from_dict(filtered_v_dict, orient="index", columns=cells.genes_name)

        df_sum_s = pd.DataFrame.from_dict(
            sum_s_dict,
            orient="index",
            columns=[
                "S0",
                "S0_norm",
                "s(t)",
                "s(t)_filtered",
                "s(t)_norm",
                "s(t)_norm_filtered",
            ],
        )

        # Store per-t results using the rounded key
        df_sum_s_dict[tk] = df_sum_s

        change_in_sum_s_per_cell_per_t = (df_sum_s["s(t)_filtered"] - df_sum_s["S0"]) / df_sum_s["S0"]
        change_in_sum_s_per_cell_all_t[tk] = change_in_sum_s_per_cell_per_t

        all_gene_level_dfs_per_t[tk] = {
            "S0": df_S0,
            "S0_norm": df_S0_norm,
            "st": df_st,
            "st_filtered": df_st_filtered,
            "st_norm": df_st_norm,
            "st_norm_filtered": df_st_norm_filtered,
            "normal_v": df_normal_v,
            "filtered_v": df_filtered_v,
        }

    with open(out_dir / "spliced_values_in_multiple_t.pkl", "wb") as f:
        pickle.dump(all_gene_level_dfs_per_t, f)

    with open(out_dir / "df_sum_s_dict.pkl", "wb") as f:
        pickle.dump(df_sum_s_dict, f)

    print(f"[stage1-d] Saved: {out_dir / 'spliced_values_in_multiple_t.pkl'}")
    print(f"[stage1-d] Saved: {out_dir / 'df_sum_s_dict.pkl'}")


    # Arrow plots
    t_fig_dir = out_dir / "t_figures"
    t_fig_dir.mkdir(parents=True, exist_ok=True)

    print("[stage1-d] Generating PCA arrow plots ...")

    for t in T2:
        tk = _tkey(t)
        S0 = all_gene_level_dfs_per_t[tk]["S0_norm"]
        St = all_gene_level_dfs_per_t[tk]["st_norm_filtered"]
        total_counts_in_cell = St.sum(axis=1)
        mean_total_counts_in_cell = total_counts_in_cell.mean()

        norm_st = St * mean_total_counts_in_cell
        scaler = StandardScaler(with_mean=True, with_std=True)
        X = scaler.fit_transform(norm_st)

        pca = PCA(n_components=3)
        pcs = pca.fit_transform(X)
        pc_St = pd.DataFrame(pcs, index=St.index, columns=["PC1", "PC2", "PC3"])

        if n_cells == "all":
            subset_cells = cells.df.copy()
        elif n_cells == "cell_type":
            subset_cells = cells.df[cells.df.cell_type == cell_type]
        else:
            subset_cells = cells.df.sample(n=int(n_cells), random_state=42)

        start_x = subset_cells["PC1"]
        start_y = subset_cells["PC2"]
        end_x = pc_St["PC1"].loc[subset_cells.index]
        end_y = pc_St["PC2"].loc[subset_cells.index]

        delta_x = end_x - start_x
        delta_y = end_y - start_y

        plt.figure(figsize=(10, 10))
        plt.scatter(start_x, start_y, c="gray", s=18, label="Original positions")
        plt.quiver(
            start_x, start_y,
            delta_x, delta_y,
            angles="xy", scale_units="xy", scale=1,
            color="blue", alpha=0.5, width=0.003,
            label="Future direction",
        )

        plt.xlabel("PC1", fontsize=32)
        plt.ylabel("PC2", fontsize=32)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.title(f"T={str(t)}", fontsize=40)

        #plt.xlim(-9, 15)
        #plt.ylim(-8, 13)
        plt.grid(True)
        ax = plt.gca()
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

        plt.tight_layout()

        
        # --- Dynamic axis limits to include full arrows ---
        all_x = np.concatenate([start_x.values, end_x.values])
        all_y = np.concatenate([start_y.values, end_y.values])

        # add 5% padding
        x_pad = 0.05 * (all_x.max() - all_x.min())
        y_pad = 0.05 * (all_y.max() - all_y.min())

        plt.xlim(all_x.min() - x_pad, all_x.max() + x_pad)
        plt.ylim(all_y.min() - y_pad, all_y.max() + y_pad)

        plt.savefig(t_fig_dir / f"t_{int(t)}-{cell_type}.tif", dpi=300)
        plt.close()

    print(f"[stage1-d] Arrow plots saved to: {t_fig_dir}")
    print("[stage1-d] Completed successfully")

