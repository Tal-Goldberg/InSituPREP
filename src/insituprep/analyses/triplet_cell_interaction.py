from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


@dataclass
class TripletParams:
    attach: float = 10.0
    alpha: float = 0.05
    n_cpus: int = 8
    refit_cooks: bool = True


def _get_coords(df: pd.DataFrame) -> np.ndarray:
    if {"X_space", "Y_space", "Z_space"}.issubset(df.columns):
        return df[["X_space", "Y_space", "Z_space"]].to_numpy(dtype=float)
    if {"X_space", "Y_space"}.issubset(df.columns):
        return df[["X_space", "Y_space"]].to_numpy(dtype=float)
    raise ValueError("Missing coordinate columns. Expected X_space,Y_space(,Z_space).")


def _find_triplet_labels(
    df: pd.DataFrame,
    attach: float,
    celltype_col: str = "cell_type",
    cell_id_col: str = "Var1",
) -> pd.DataFrame:
    """
    For each cell, find neighbor cell types within 'attach' radius,
    keep only cases with up to 3 unique types total (<=2 commas).
    Returns: columns [cell, type, friends]
    """
    coords = _get_coords(df)
    tree = cKDTree(coords)

    cell_ids = df[cell_id_col].astype(str).to_numpy()
    types = df[celltype_col].astype(str).to_numpy()

    friends_list: List[str] = []
    keep_cell: List[str] = []
    keep_type: List[str] = []

    for i in range(len(df)):
        idxs = tree.query_ball_point(coords[i], r=attach)
        if len(idxs) <= 1:
            continue

        uniq_types = np.unique(types[idxs])
        # remove "nan" strings if exist
        uniq_types = uniq_types[uniq_types != "nan"]

        if len(uniq_types) == 0:
            continue

        friends = ",".join(sorted(uniq_types))
        # Keep only up to triplets (<= 3 types => <=2 commas)
        if friends.count(",") <= 2:
            keep_cell.append(cell_ids[i])
            keep_type.append(types[i])
            friends_list.append(friends)

    out = pd.DataFrame({"cell": keep_cell, "type": keep_type, "friends": friends_list})
    out = out[out["type"].astype(str) != "nan"].copy()
    out = out[out["friends"].astype(str) != "nan"].copy()
    return out


def _build_deseq_inputs(
    df: pd.DataFrame,
    triplets_df: pd.DataFrame,
    genes: List[str],
    cell_id_col: str = "Var1",
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Reproduce your A/B group building logic in a safer/clearer way.
    NOTE: This logic is kept close to your script; you can refine later.
    """
    # +1 pseudocount as you did
    m = df.copy()
    for g in genes:
        if g not in m.columns:
            raise ValueError(f"Gene '{g}' not found in labels/counts table.")
    m[genes] = m[genes].to_numpy() + 1

    counts: Dict[str, pd.DataFrame] = {}
    metadata: Dict[str, pd.DataFrame] = {}

    # Work per focal type i
    for i_type in sorted(triplets_df["type"].unique()):
        sub = triplets_df[triplets_df["type"] == i_type].copy()

        # only friends with at least 2 types (has comma)
        candidate_friends = sorted(set([s for s in sub["friends"].unique() if "," in s]))
        for k in candidate_friends:
            A = sub[sub["friends"] == k]
            if A.empty:
                continue

            # In your original code you tried to create B via k+','+j
            # This is order-sensitive. We'll keep the same idea but check both orders.
            for j_type in sorted(triplets_df["type"].unique()):
                if j_type == i_type:
                    continue

                k1 = f"{k},{j_type}"
                k2 = f"{j_type},{k}"
                B = sub[(sub["friends"] == k1) | (sub["friends"] == k2)]
                if B.empty:
                    continue

                key = f"{i_type}-{k}-{j_type}"

                cells_A = A["cell"].astype(str).tolist()
                cells_B = B["cell"].astype(str).tolist()

                X = pd.concat(
                    [
                        m[m[cell_id_col].astype(str).isin(cells_A)][genes],
                        m[m[cell_id_col].astype(str).isin(cells_B)][genes],
                    ],
                    axis=0,
                )

                X.index = [f"Sample{idx}" for idx in range(X.shape[0])]

                y = pd.DataFrame({"condition": ["G"] * len(cells_A) + ["A"] * len(cells_B)})
                y.index = X.index

                counts[key] = X
                metadata[key] = y

    return counts, metadata


def _run_deseq_for_group(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    params: TripletParams,
) -> pd.DataFrame:
    inference = DefaultInference(n_cpus=params.n_cpus)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=meta_df,
        design="~condition",
        refit_cooks=params.refit_cooks,
        inference=inference,
    )
    dds.deseq2()

    ds = DeseqStats(dds, contrast=["condition", "G", "A"], inference=inference)
    ds.summary()
    res = ds.results_df.copy()

    # filter
    res = res[np.isfinite(res["pvalue"])]
    res = res[res["pvalue"] < params.alpha]
    return res


def run_triplet_cell_interaction(
    labels_df: pd.DataFrame,
    genes: List[str],
    params: TripletParams,
    tissue: Optional[str] = None,
    tissue_col: str = "tissue",
    celltype_col: str = "cell_type",
    cell_id_col: str = "Var1",
) -> pd.DataFrame:
    """
    Main entry point:
    - if tissue is provided and tissue_col exists => subset
    - compute triplet labels
    - build DESeq2 inputs and run per group
    """
    df = labels_df.copy()
    df[cell_id_col] = df[cell_id_col].astype(str)
    if tissue is not None:
        if tissue_col not in df.columns:
            raise ValueError(f"--tissue was provided but column '{tissue_col}' not found.")
        df = df[df[tissue_col].astype(str) == str(tissue)].copy()

    trip = _find_triplet_labels(df, attach=params.attach, celltype_col=celltype_col, cell_id_col=cell_id_col)
    if trip.empty:
        return pd.DataFrame()

    counts_map, meta_map = _build_deseq_inputs(df, trip, genes=genes, cell_id_col=cell_id_col)
    if not counts_map:
        return pd.DataFrame()

    out_rows = []
    for key, X in counts_map.items():
        y = meta_map[key]
        try:
            res = _run_deseq_for_group(X, y, params=params)
            if res is None or res.empty:
                continue
            res = res.copy()
            res["triplet_group"] = key
            out_rows.append(res)
        except Exception:
            # keep behavior similar to your "except: pass" but without killing the run
            continue

    if not out_rows:
        return pd.DataFrame()

    result = pd.concat(out_rows, axis=0)

    # Split 'triplet_group' into 3 columns: focal_type, pair_types, neighbor_type
    # Use maxsplit=2 so middle part may contain commas (e.g. "A,B").
    parts = result["triplet_group"].astype(str).str.split("-", n=2, expand=True)
    parts.columns = ["triplet_focal_type", "triplet_pair_types", "triplet_neighbor_type"]

    # join back to result (keeps original 'triplet_group' as well)
    result = pd.concat([result.reset_index(drop=True), parts.reset_index(drop=True)], axis=1)

    # Optional: if you prefer to drop the original column, uncomment:
    result = result.drop(columns=["triplet_group"])

    return result


    return pd.concat(out_rows, axis=0)
