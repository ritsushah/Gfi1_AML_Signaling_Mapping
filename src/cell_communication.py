"""Stage 3 – Ligand–receptor interaction scoring between leukemic and niche cells."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import LR_PAIRS, SENDER_TYPES, RECEIVER_TYPES, TABLES


def score_interactions(
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    lr_pairs: list[tuple[str, str]] | None = None,
    min_expr: float = 0.05,
) -> pd.DataFrame:
    """
    Simple but interpretable mean-expression product score for each
    (ligand, receptor, sender, receiver) combination.

    For production, replace with CellPhoneDB / LIANA / NicheNet.
    """
    if lr_pairs is None:
        lr_pairs = LR_PAIRS

    records = []
    cell_types = meta["cell_type"].unique()

    for ligand, receptor in lr_pairs:
        if ligand not in expr.columns or receptor not in expr.columns:
            continue

        for sender in SENDER_TYPES:
            if sender not in cell_types:
                continue
            sender_mask = meta["cell_type"] == sender
            # Optionally further restrict senders to GFI1-Low
            if "GFI1_status" in meta.columns:
                sender_mask = sender_mask & (meta["GFI1_status"] == "Low")

            if sender_mask.sum() < 5:
                continue
            lig_mean = expr.loc[sender_mask, ligand].mean()

            for receiver in RECEIVER_TYPES:
                if receiver not in cell_types:
                    continue
                receiver_mask = meta["cell_type"] == receiver
                if receiver_mask.sum() < 5:
                    continue
                rec_mean = expr.loc[receiver_mask, receptor].mean()

                if lig_mean >= min_expr and rec_mean >= min_expr:
                    score = lig_mean * rec_mean
                    records.append({
                        "ligand": ligand,
                        "receptor": receptor,
                        "sender": sender,
                        "receiver": receiver,
                        "ligand_mean": round(lig_mean, 4),
                        "receptor_mean": round(rec_mean, 4),
                        "interaction_score": round(score, 4),
                    })

    df = pd.DataFrame(records)
    if len(df) > 0:
        df = df.sort_values("interaction_score", ascending=False).reset_index(drop=True)
    return df


def save_communication_results(interactions: pd.DataFrame) -> None:
    path = TABLES / "lr_interaction_scores.csv"
    interactions.to_csv(path, index=False)
    print(f"[Comm] Interaction table → {path}  ({len(interactions)} pairs)")
