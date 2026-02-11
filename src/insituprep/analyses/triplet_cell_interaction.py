from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import cdist

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


@dataclass
class TripletParams:
    attach: float = 10.0
    n_cpus: int = 8
    refit_cooks: bool = True


def _get_coords(df: pd.DataFrame) -> np.ndarray:
   
    if {"X_space", "Y_space", "Z_space"}.issubset(df.columns):
        return df[["X_space", "Y_space", "Z_space"]].to_numpy(dtype=float)
    if {"X_space", "Y_space"}.issubset(df.columns):
        return df[["X_space", "Y_space"]].to_numpy(dtype=float)
    raise ValueError("Missing coordinate columns. Expected X_space,Y_space(,Z_space).")


def _find_triplet_labels_cdist(
    df: pd.DataFrame,
    attach: float,
    celltype_col: str = "cell_type",
    cell_id_col: str = "Var1",
) -> pd.DataFrame:
    """
    Implements the UPDATED logic (new code):
    - For each cell i: compute distances d = cdist([z[i]], z)
    - If count(d < attach) > 1: keep the cell
    - friends := ','.join(sorted(unique(types within attach)))
    - Filter out rows where type is 'nan'
    - Keep only friends strings with <= 2 commas (<= 3 unique types)
    """
    z = _get_coords(df)
    cell_ids = df[cell_id_col].astype(str).to_numpy()
    types = df[celltype_col].astype(str).to_numpy()

    kept_cell: List[str] = []
    kept_type: List[str] = []
    kept_friends: List[str] = []

    for i in range(df.shape[0]):
        d = cdist(np.array([z[i]]), z)  # shape (1, n)
        mask = (d < attach).ravel()

        # בדיוק כמו בקוד החדש: תנאי "יותר מאחד" (כולל self)
        if np.sum(mask) > 1:
            kept_cell.append(cell_ids[i])
            kept_type.append(types[i])

            friends_types = np.unique(types[mask])
            friends_types = friends_types[friends_types != "nan"]
            kept_friends.append(",".join(np.sort(friends_types)))

    c = pd.DataFrame({"cell": kept_cell, "type": kept_type, "friends": kept_friends})

    c = c[c["type"].astype(str) != "nan"].copy()

    # <=2 commas (triplets)
    comma_count = c["friends"].str.len() - c["friends"].str.replace(",", "", regex=False).str.len()
    c = c[comma_count <= 2].copy()

    return c


def _build_deseq_inputs_updated(
    labels_df: pd.DataFrame,
    triplets_df: pd.DataFrame,
    genes: List[str],
    cell_id_col: str = "Var1",
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Implements the UPDATED grouping logic (new code), including:
    - m = labels_df ; m[genes] += 1  (pseudocount)
    - For each primary type i in unique(triplets_df.type):
        x = triplets_df[triplets_df.type == i]
        for each k in unique(x.friends) where k contains comma:
            A = x[x.friends == k]
            for each j in unique(triplets_df.type):
                B = x[x.friends == k + ',' + j]   (ORDER-SENSITIVE, no reverse!)
                if A and B non-empty: create counts/metadata with key f"{i}-{k}-{j}"
    """
    m = labels_df.copy()

    # validate genes exist
    missing = [g for g in genes if g not in m.columns]
    if missing:
        raise ValueError(f"Missing gene columns in labels table: {missing[:10]}{'...' if len(missing) > 10 else ''}")

    # exactly like new code: in-place add 1
    m.loc[:, genes] = m[genes].to_numpy() + 1

    counts: Dict[str, pd.DataFrame] = {}
    metadata: Dict[str, pd.DataFrame] = {}

    all_types = np.unique(triplets_df["type"].astype(str).to_numpy())

    for i_type in np.unique(triplets_df["type"].astype(str).to_numpy()):
        x = triplets_df[triplets_df["type"].astype(str) == str(i_type)].copy()

        # only k with at least one comma
        ks = np.unique(x["friends"][x["friends"].astype(str).str.find(",") > -1].astype(str).to_numpy())
        for k in ks:
            A = x[x["friends"].astype(str) == str(k)]
            if A.shape[0] == 0:
                continue

            for j_type in all_types:
                # IMPORTANT: exactly like new code (order-sensitive)
                B = x[x["friends"].astype(str) == f"{k},{j_type}"]
                if A.shape[0] > 0 and B.shape[0] > 0:
                    key = f"{i_type}-{k}-{j_type}"

                    s_counts = pd.concat(
                        [
                            m[m[cell_id_col].astype(str).isin(A["cell"].astype(str))][genes],
                            m[m[cell_id_col].astype(str).isin(B["cell"].astype(str))][genes],
                        ],
                        axis=0,
                    )
                    s_counts.index = "Sample" + pd.Series(np.arange(s_counts.shape[0])).astype(str).to_numpy()
                    counts[key] = s_counts

                    s_meta = pd.DataFrame({"condition": ["G"] * A.shape[0] + ["A"] * B.shape[0]})
                    s_meta.index = "Sample" + pd.Series(np.arange(s_meta.shape[0])).astype(str).to_numpy()
                    metadata[key] = s_meta

    return counts, metadata


def _run_deseq_unfiltered(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    params: TripletParams,
) -> pd.DataFrame:
    """
    IMPORTANT: matches UPDATED behavior:
    - run DESeq2
    - return ds.results_df WITHOUT filtering by alpha/pvalue here
    """
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
    return ds.results_df.copy()


def _annotate_results_from_key(res: pd.DataFrame, key: str) -> pd.DataFrame:
    """
    Matches UPDATED annotation logic exactly:
    Z['primary type'] = split[0]
    Z['neighbor type'] = split[1].replace(split[0],'').replace(',','')
    Z['neighbor added type'] = split[2]
    """
    parts = str(key).split("-")
    primary = parts[0] if len(parts) > 0 else ""
    pair = parts[1] if len(parts) > 1 else ""
    added = parts[2] if len(parts) > 2 else ""

    out = res.copy()
    out["primary type"] = primary
    out["neighbor type"] = str(pair).replace(primary, "").replace(",", "")
    out["neighbor added type"] = added
    return out


def run_triplet_cell_interaction(
    labels_df: pd.DataFrame,
    genes: List[str],
    params: TripletParams,
    pv_threshold: float,
    tissue: Optional[str] = None,
    tissue_col: str = "tissue",
    celltype_col: str = "cell_type",
    cell_id_col: str = "Var1",
) -> pd.DataFrame:
    """
    Package-style entry point implementing UPDATED script behavior:

    If tissue_col NOT in labels_df:
        - run once
        - DO NOT add padj / DO NOT filter (like script's "no tissue column" branch)
        - return concatenated results

    If tissue_col exists:
        - If tissue is provided: run only that tissue
        - Else: run per tissue and concatenate
        - For EACH tissue: set non-finite pvalue -> 1, compute padj via scipy.stats.false_discovery_control,
          filter padj < pv_threshold (exactly like script)
    """
    df = labels_df.copy()
    df[cell_id_col] = df[cell_id_col].astype(str)

    def _run_one(sub_df: pd.DataFrame) -> pd.DataFrame:
        trip = _find_triplet_labels_cdist(
            sub_df,
            attach=params.attach,
            celltype_col=celltype_col,
            cell_id_col=cell_id_col,
        )
        if trip.empty:
            return pd.DataFrame()

        counts_map, meta_map = _build_deseq_inputs_updated(
            labels_df=sub_df,
            triplets_df=trip,
            genes=genes,
            cell_id_col=cell_id_col,
        )
        if not counts_map:
            return pd.DataFrame()

        all_res: List[pd.DataFrame] = []
        for key, X in counts_map.items():
            y = meta_map[key]
            try:
                Z = _run_deseq_unfiltered(X, y, params=params)
                Z = _annotate_results_from_key(Z, key)
                all_res.append(Z)
            except Exception:
                # keep "except: pass" behavior from script
                continue

        if not all_res:
            return pd.DataFrame()

        return pd.concat(all_res, axis=0)

    # === Case 1: no tissue column ===
    if tissue_col not in df.columns:
        return _run_one(df)

    # === Case 2: tissue column exists ===
    if tissue is not None:
        sub = df[df[tissue_col].astype(str) == str(tissue)].copy()
        sub.index = np.arange(sub.shape[0])
        z = _run_one(sub)
        if z.empty:
            return z

        # UPDATED script behavior: pvalue nonfinite -> 1; padj via false_discovery_control; filter padj < pv_threshold
        z["pvalue"] = np.where(np.isfinite(z["pvalue"]), z["pvalue"], 1)
        z["padj"] = scipy.stats.false_discovery_control(z["pvalue"])
        z = z[z["padj"] < float(pv_threshold)].copy()
        z[tissue_col] = str(tissue)
        return z

    # run all tissues and concat (each tissue filtered independently, like script files-per-tissue)
    out: List[pd.DataFrame] = []
    for tis in np.unique(df[tissue_col].astype(str).to_numpy()):
        sub = df[df[tissue_col].astype(str) == str(tis)].copy()
        sub.index = np.arange(sub.shape[0])

        z = _run_one(sub)
        if z.empty:
            continue

        z["pvalue"] = np.where(np.isfinite(z["pvalue"]), z["pvalue"], 1)
        z["padj"] = scipy.stats.false_discovery_control(z["pvalue"])
        z = z[z["padj"] < float(pv_threshold)].copy()
        if z.empty:
            continue

        z[tissue_col] = str(tis)
        out.append(z)

    if not out:
        return pd.DataFrame()

    return pd.concat(out, axis=0)
