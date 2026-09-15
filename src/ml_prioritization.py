"""Stage 2 – Machine-learning prioritization of secreted candidates using clinical outcome."""

from __future__ import annotations

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import xgboost as xgb

from .config import RANDOM_SEED, TABLES, MODELS


def train_risk_model(
    expr: pd.DataFrame,
    clinical: pd.DataFrame,
    label_col: str = "high_risk",
) -> tuple[Pipeline, pd.Series, dict]:
    """
    Train an XGBoost classifier to predict high-risk status from gene expression.
    Returns fitted pipeline, feature-importance Series, and performance metrics.
    """
    # Align samples
    common = expr.index.intersection(clinical.index)
    X = expr.loc[common]
    y = clinical.loc[common, label_col]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=RANDOM_SEED, stratify=y
    )

    model = xgb.XGBClassifier(
        n_estimators=300,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.85,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=RANDOM_SEED,
        eval_metric="logloss",
        n_jobs=-1,
    )

    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", model),
    ])
    pipe.fit(X_train, y_train)

    # Metrics
    proba = pipe.predict_proba(X_test)[:, 1]
    metrics = {
        "test_auc": float(roc_auc_score(y_test, proba)),
        "test_ap": float(average_precision_score(y_test, proba)),
        "n_train": len(X_train),
        "n_test": len(X_test),
    }

    # Cross-validated AUC on full data
    cv_scores = cross_val_score(pipe, X, y, cv=5, scoring="roc_auc", n_jobs=-1)
    metrics["cv_auc_mean"] = float(cv_scores.mean())
    metrics["cv_auc_std"] = float(cv_scores.std())

    # Feature importance (from the underlying booster)
    clf = pipe.named_steps["clf"]
    importance = pd.Series(
        clf.feature_importances_, index=X.columns, name="importance"
    ).sort_values(ascending=False)

    return pipe, importance, metrics


def save_ml_results(pipe, importance: pd.Series, metrics: dict) -> None:
    importance.to_csv(TABLES / "ml_feature_importance.csv")
    joblib.dump(pipe, MODELS / "xgboost_risk_model.joblib")

    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(TABLES / "ml_performance.csv", index=False)

    print(f"[ML] Feature importance → {TABLES / 'ml_feature_importance.csv'}")
    print(f"[ML] Model saved        → {MODELS / 'xgboost_risk_model.joblib'}")
    print(f"[ML] Test AUC = {metrics['test_auc']:.3f} | CV AUC = {metrics['cv_auc_mean']:.3f} ± {metrics['cv_auc_std']:.3f}")
