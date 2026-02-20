"""QSAR model building and prediction tool."""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.qsar_model")


class QSARModelTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="qsar_model",
            description=(
                "Build a QSAR (Quantitative Structure-Activity Relationship) model from "
                "molecular fingerprints and activity data. Supports Random Forest and XGBoost "
                "regressors. Returns cross-validation metrics (R-squared, RMSE, MAE) and "
                "optionally predicts activities for new compounds."
            ),
            parameters=[
                ToolParameter(
                    name="smiles_list",
                    type=ParameterType.ARRAY,
                    description="List of SMILES strings for training compounds",
                    items_type=ParameterType.STRING,
                ),
                ToolParameter(
                    name="activities",
                    type=ParameterType.ARRAY,
                    description="List of activity values corresponding to each SMILES",
                    items_type=ParameterType.NUMBER,
                ),
                ToolParameter(
                    name="model_type",
                    type=ParameterType.STRING,
                    description="Machine learning model to use",
                    required=False,
                    default="random_forest",
                    enum=["random_forest", "xgboost"],
                ),
                ToolParameter(
                    name="fingerprint_type",
                    type=ParameterType.STRING,
                    description="Type of molecular fingerprint for feature generation",
                    required=False,
                    default="ecfp4",
                    enum=["ecfp4", "maccs"],
                ),
                ToolParameter(
                    name="predict_smiles",
                    type=ParameterType.ARRAY,
                    description="Optional list of SMILES strings to predict activity for",
                    required=False,
                    items_type=ParameterType.STRING,
                ),
            ],
            category="drugdiscovery",
            is_long_running=True,
        )

    def check_available(self) -> bool:
        try:
            from rdkit import Chem  # noqa: F401
            from sklearn.ensemble import RandomForestRegressor  # noqa: F401

            return True
        except ImportError:
            return False

    def _smiles_to_fingerprint(self, smiles: str, fp_type: str) -> np.ndarray | None:
        """Convert a SMILES string to a fingerprint bit vector."""
        from rdkit import Chem
        from rdkit.Chem import AllChem, MACCSkeys

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        if fp_type == "ecfp4":
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        else:
            fp = MACCSkeys.GenMACCSKeys(mol)

        arr = np.zeros(fp.GetNumBits(), dtype=np.int8)
        for bit in fp.GetOnBits():
            arr[bit] = 1
        return arr

    def _build_fingerprint_matrix(
        self, smiles_list: list[str], fp_type: str
    ) -> tuple[np.ndarray, list[int], list[str]]:
        """Build fingerprint matrix, returning valid indices and warnings."""
        fps = []
        valid_indices = []
        warnings: list[str] = []

        for i, smi in enumerate(smiles_list):
            fp = self._smiles_to_fingerprint(smi, fp_type)
            if fp is None:
                warnings.append(f"Could not parse SMILES at index {i}: {smi}")
            else:
                fps.append(fp)
                valid_indices.append(i)

        X = np.array(fps) if fps else np.empty((0, 0))
        return X, valid_indices, warnings

    async def execute(self, **kwargs) -> ToolResult:
        smiles_list: list[str] = kwargs["smiles_list"]
        activities: list[float] = kwargs["activities"]
        model_type: str = kwargs.get("model_type", "random_forest")
        fingerprint_type: str = kwargs.get("fingerprint_type", "ecfp4")
        predict_smiles: list[str] | None = kwargs.get("predict_smiles")

        # Validate inputs
        if len(smiles_list) != len(activities):
            return ToolResult(
                success=False,
                error=f"Length mismatch: {len(smiles_list)} SMILES vs {len(activities)} activities",
                summary="SMILES list and activities list must have the same length.",
            )

        if len(smiles_list) < 5:
            return ToolResult(
                success=False,
                error="At least 5 compounds are needed for cross-validation",
                summary="Insufficient data: need at least 5 compounds for QSAR modeling.",
            )

        try:
            from rdkit import Chem  # noqa: F401
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.model_selection import cross_val_score
        except ImportError:
            return ToolResult(
                success=False,
                error="Required packages not installed. Install rdkit and scikit-learn.",
                summary="Missing dependencies: rdkit and/or scikit-learn",
            )

        # Build fingerprint matrix
        X, valid_indices, warnings = self._build_fingerprint_matrix(smiles_list, fingerprint_type)
        y = np.array([activities[i] for i in valid_indices])

        if len(valid_indices) < 5:
            return ToolResult(
                success=False,
                error=f"Only {len(valid_indices)} valid molecules after parsing (need >= 5)",
                summary="Too few valid molecules for QSAR modeling after SMILES parsing.",
                warnings=warnings,
            )

        # Build model
        if model_type == "xgboost":
            try:
                from xgboost import XGBRegressor

                model = XGBRegressor(
                    n_estimators=100,
                    max_depth=6,
                    learning_rate=0.1,
                    random_state=42,
                    verbosity=0,
                )
            except ImportError:
                return ToolResult(
                    success=False,
                    error="xgboost is not installed. Install with: pip install xgboost",
                    summary="XGBoost not available",
                )
        else:
            model = RandomForestRegressor(
                n_estimators=100,
                max_depth=None,
                random_state=42,
                n_jobs=-1,
            )

        # Cross-validation
        n_folds = min(5, len(valid_indices))

        r2_scores = cross_val_score(model, X, y, cv=n_folds, scoring="r2")
        neg_mse_scores = cross_val_score(model, X, y, cv=n_folds, scoring="neg_mean_squared_error")
        neg_mae_scores = cross_val_score(model, X, y, cv=n_folds, scoring="neg_mean_absolute_error")

        rmse_scores = np.sqrt(-neg_mse_scores)
        mae_scores = -neg_mae_scores

        cv_metrics = {
            "r2_mean": round(float(r2_scores.mean()), 4),
            "r2_std": round(float(r2_scores.std()), 4),
            "rmse_mean": round(float(rmse_scores.mean()), 4),
            "rmse_std": round(float(rmse_scores.std()), 4),
            "mae_mean": round(float(mae_scores.mean()), 4),
            "mae_std": round(float(mae_scores.std()), 4),
            "n_folds": n_folds,
        }

        # Fit final model on all data
        model.fit(X, y)

        # Feature importance (top 20 bits)
        importances = model.feature_importances_
        top_indices = np.argsort(importances)[-20:][::-1]
        feature_importance = [
            {"bit": int(idx), "importance": round(float(importances[idx]), 6)}
            for idx in top_indices
            if importances[idx] > 0
        ]

        # Predictions for new compounds
        predictions = None
        pred_warnings: list[str] = []
        if predict_smiles:
            pred_X, pred_valid, pred_warns = self._build_fingerprint_matrix(
                predict_smiles, fingerprint_type
            )
            pred_warnings = pred_warns

            if len(pred_valid) > 0:
                pred_values = model.predict(pred_X)
                predictions = []
                pred_idx = 0
                for i, smi in enumerate(predict_smiles):
                    if i in [pred_valid[j] for j in range(len(pred_valid))]:
                        predictions.append({
                            "smiles": smi,
                            "predicted_activity": round(float(pred_values[pred_idx]), 4),
                            "index": i,
                        })
                        pred_idx += 1
                    else:
                        predictions.append({
                            "smiles": smi,
                            "predicted_activity": None,
                            "index": i,
                            "error": "Could not parse SMILES",
                        })

        all_warnings = warnings + pred_warnings

        data: dict[str, Any] = {
            "model_type": model_type,
            "fingerprint_type": fingerprint_type,
            "training_size": len(valid_indices),
            "cv_metrics": cv_metrics,
            "feature_importance": feature_importance,
        }
        if predictions is not None:
            data["predictions"] = predictions

        # Build summary
        r2 = cv_metrics["r2_mean"]
        rmse = cv_metrics["rmse_mean"]
        quality = "good" if r2 > 0.7 else "moderate" if r2 > 0.4 else "poor"

        summary_lines = [
            f"QSAR model ({model_type}, {fingerprint_type.upper()}): "
            f"trained on {len(valid_indices)} compounds.",
            f"Cross-validation ({n_folds}-fold): R2={r2:.4f} (+/-{cv_metrics['r2_std']:.4f}), "
            f"RMSE={rmse:.4f}, MAE={cv_metrics['mae_mean']:.4f}. Model quality: {quality}.",
        ]
        if predictions:
            valid_preds = [p for p in predictions if p["predicted_activity"] is not None]
            summary_lines.append(f"Predictions made for {len(valid_preds)} compound(s).")
            for p in valid_preds[:5]:
                summary_lines.append(f"  - {p['smiles']}: {p['predicted_activity']}")

        if r2 < 0.4:
            all_warnings.append("Model has poor predictive performance (R2 < 0.4). Use with caution.")

        return ToolResult(
            success=True,
            data=data,
            summary="\n".join(summary_lines),
            warnings=all_warnings,
        )
