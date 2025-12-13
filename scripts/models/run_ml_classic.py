#!/usr/bin/env python3
"""
Platelet Classification - Classic Machine Learning Pipeline
============================================================
Predicts patient survival (S vs NS) using single-cell platelet gene expression.

Models: XGBoost, Random Forest, Gradient Boosting, Logistic Regression, SVM

Usage:
    cd /bigdata/godziklab/shared/Xinru/302005
    python3 PlateletSubpop-ML-ScTranscriptomics/scripts/models/run_ml_classic.py

Output files saved to current directory:
    - model_comparison_summary.csv
    - roc_curves_comparison.png
    - confusion_matrices.png
    - feature_importance_xgboost.csv/png
    - feature_importance_rf.csv/png
    - consensus_top_genes.csv
    - cv_results_boxplot.png
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for HPCC
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (accuracy_score, roc_curve, roc_auc_score,
                             confusion_matrix, classification_report)
from xgboost import XGBClassifier
import warnings
import os

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
CONFIG = {
    'data_path': '302005_platelet_harmony_integrated.h5ad',
    'target_col': 'Death',
    'target_classes': ['S', 'NS'],  # Survivor vs Non-Survivor
    'n_top_genes': 500,
    'test_size': 0.2,
    'cv_folds': 5,
    'random_state': 42
}

# ============================================================
# Main Pipeline
# ============================================================
def main():
    print("=" * 60)
    print("PLATELET CLASSIFICATION - CLASSIC ML PIPELINE")
    print("=" * 60)

    # ----------------------------------------------------------
    # 1. Load Data
    # ----------------------------------------------------------
    print("\n[1/7] Loading data...")

    if not os.path.exists(CONFIG['data_path']):
        print(f"ERROR: Data file not found: {CONFIG['data_path']}")
        print("Make sure you run this from the directory containing the h5ad file")
        return

    adata = sc.read_h5ad(CONFIG['data_path'])
    print(f"  Original: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # Filter to target classes
    mask = adata.obs[CONFIG['target_col']].isin(CONFIG['target_classes'])
    adata = adata[mask].copy()
    print(f"  Filtered ({' vs '.join(CONFIG['target_classes'])}): {adata.n_obs:,} cells")

    class_counts = adata.obs[CONFIG['target_col']].value_counts()
    for cls, count in class_counts.items():
        print(f"    - {cls}: {count:,} ({count/len(adata)*100:.1f}%)")

    # ----------------------------------------------------------
    # 2. Prepare Features
    # ----------------------------------------------------------
    print("\n[2/7] Preparing features...")

    # Extract expression matrix
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(adata.obs[CONFIG['target_col']])
    print(f"  Classes: {dict(zip(le.classes_, range(len(le.classes_))))}")

    # Feature selection: top variable genes
    gene_vars = np.var(X, axis=0)
    top_idx = np.argsort(gene_vars)[-CONFIG['n_top_genes']:]
    X_selected = X[:, top_idx]
    gene_names = adata.var_names[top_idx].tolist()
    print(f"  Selected top {CONFIG['n_top_genes']} variable genes")
    print(f"  Feature matrix: {X_selected.shape}")

    # ----------------------------------------------------------
    # 3. Train/Test Split
    # ----------------------------------------------------------
    print("\n[3/7] Splitting data...")

    X_train, X_test, y_train, y_test = train_test_split(
        X_selected, y,
        test_size=CONFIG['test_size'],
        random_state=CONFIG['random_state'],
        stratify=y
    )

    # Scale for SVM and Logistic Regression
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    print(f"  Train: {len(y_train):,} samples")
    print(f"  Test:  {len(y_test):,} samples")

    # ----------------------------------------------------------
    # 4. Define Models
    # ----------------------------------------------------------
    print("\n[4/7] Defining models...")

    models = {
        'XGBoost': XGBClassifier(
            n_estimators=100, max_depth=6, learning_rate=0.1,
            random_state=CONFIG['random_state'],
            use_label_encoder=False, eval_metric='logloss',
            n_jobs=-1
        ),
        'Random Forest': RandomForestClassifier(
            n_estimators=200, max_depth=10,
            random_state=CONFIG['random_state'], n_jobs=-1
        ),
        'Gradient Boosting': GradientBoostingClassifier(
            n_estimators=100, max_depth=5, learning_rate=0.1,
            random_state=CONFIG['random_state']
        ),
        'Logistic Regression': LogisticRegression(
            max_iter=1000, random_state=CONFIG['random_state'], n_jobs=-1
        ),
        'SVM': SVC(
            kernel='rbf', probability=True,
            random_state=CONFIG['random_state']
        )
    }
    print(f"  Models: {list(models.keys())}")

    # ----------------------------------------------------------
    # 5. Train and Evaluate
    # ----------------------------------------------------------
    print("\n[5/7] Training and evaluating models...")

    results = {}
    cv = StratifiedKFold(n_splits=CONFIG['cv_folds'], shuffle=True,
                         random_state=CONFIG['random_state'])

    for name, model in models.items():
        print(f"\n  >> {name}")

        # Select appropriate data (scaled for SVM/LR)
        if name in ['SVM', 'Logistic Regression']:
            X_tr, X_te = X_train_scaled, X_test_scaled
        else:
            X_tr, X_te = X_train, X_test

        # Train
        model.fit(X_tr, y_train)

        # Predict
        y_pred = model.predict(X_te)
        y_prob = model.predict_proba(X_te)[:, 1]

        # Cross-validation
        cv_scores = cross_val_score(model, X_tr, y_train, cv=cv, scoring='roc_auc')

        # Metrics
        acc = accuracy_score(y_test, y_pred)
        auc_score = roc_auc_score(y_test, y_prob)

        results[name] = {
            'model': model,
            'y_pred': y_pred,
            'y_prob': y_prob,
            'accuracy': acc,
            'auc': auc_score,
            'cv_scores': cv_scores
        }

        print(f"     Accuracy: {acc:.4f}")
        print(f"     AUC:      {auc_score:.4f}")
        print(f"     CV-AUC:   {cv_scores.mean():.4f} (+/- {cv_scores.std()*2:.4f})")

    # ----------------------------------------------------------
    # 6. Generate Results Summary
    # ----------------------------------------------------------
    print("\n[6/7] Generating summary...")

    summary_data = []
    for name, res in results.items():
        summary_data.append({
            'Model': name,
            'Accuracy': round(res['accuracy'], 4),
            'AUC': round(res['auc'], 4),
            'CV_AUC_Mean': round(res['cv_scores'].mean(), 4),
            'CV_AUC_Std': round(res['cv_scores'].std(), 4)
        })

    summary_df = pd.DataFrame(summary_data).sort_values('AUC', ascending=False)
    summary_df.to_csv('model_comparison_summary.csv', index=False)

    print("\n" + "=" * 60)
    print("MODEL PERFORMANCE SUMMARY")
    print("=" * 60)
    print(summary_df.to_string(index=False))

    # ----------------------------------------------------------
    # 7. Generate Visualizations
    # ----------------------------------------------------------
    print("\n[7/7] Generating visualizations...")

    # Color palette
    colors = {
        'XGBoost': '#E41A1C',
        'Random Forest': '#4DAF4A',
        'Gradient Boosting': '#377EB8',
        'Logistic Regression': '#FF7F00',
        'SVM': '#984EA3'
    }

    # --- ROC Curves ---
    plt.figure(figsize=(10, 8))
    for name, res in results.items():
        fpr, tpr, _ = roc_curve(y_test, res['y_prob'])
        plt.plot(fpr, tpr, color=colors[name], lw=2,
                 label=f"{name} (AUC = {res['auc']:.3f})")

    plt.plot([0, 1], [0, 1], 'k--', lw=1, label='Random (AUC = 0.5)')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curves - Survival Prediction (S vs NS)', fontsize=14)
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('roc_curves_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: roc_curves_comparison.png")

    # --- Confusion Matrices ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, (name, res) in enumerate(results.items()):
        cm = confusion_matrix(y_test, res['y_pred'])
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[idx],
                    xticklabels=le.classes_, yticklabels=le.classes_,
                    annot_kws={'size': 14})
        axes[idx].set_title(f"{name}\nAcc: {res['accuracy']:.3f}, AUC: {res['auc']:.3f}",
                           fontsize=11)
        axes[idx].set_xlabel('Predicted', fontsize=10)
        axes[idx].set_ylabel('Actual', fontsize=10)

    axes[-1].axis('off')
    plt.tight_layout()
    plt.savefig('confusion_matrices.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: confusion_matrices.png")

    # --- CV Box Plot ---
    cv_data = {name: res['cv_scores'] for name, res in results.items()}
    cv_df = pd.DataFrame(cv_data)

    plt.figure(figsize=(10, 6))
    cv_df.boxplot(grid=False)
    plt.ylabel('AUC Score', fontsize=12)
    plt.title(f'{CONFIG["cv_folds"]}-Fold Cross-Validation Results', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('cv_results_boxplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: cv_results_boxplot.png")

    # --- Feature Importance (XGBoost) ---
    xgb_imp = pd.DataFrame({
        'Gene': gene_names,
        'Importance': results['XGBoost']['model'].feature_importances_
    }).sort_values('Importance', ascending=False)
    xgb_imp.to_csv('feature_importance_xgboost.csv', index=False)

    plt.figure(figsize=(10, 8))
    top20 = xgb_imp.head(20)
    plt.barh(range(len(top20)), top20['Importance'].values, color='steelblue')
    plt.yticks(range(len(top20)), top20['Gene'].values)
    plt.xlabel('Importance Score', fontsize=12)
    plt.title('XGBoost Feature Importance (Top 20 Genes)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig('feature_importance_xgboost.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: feature_importance_xgboost.png")

    # --- Feature Importance (Random Forest) ---
    rf_imp = pd.DataFrame({
        'Gene': gene_names,
        'Importance': results['Random Forest']['model'].feature_importances_
    }).sort_values('Importance', ascending=False)
    rf_imp.to_csv('feature_importance_rf.csv', index=False)

    plt.figure(figsize=(10, 8))
    top20_rf = rf_imp.head(20)
    plt.barh(range(len(top20_rf)), top20_rf['Importance'].values, color='forestgreen')
    plt.yticks(range(len(top20_rf)), top20_rf['Gene'].values)
    plt.xlabel('Importance Score', fontsize=12)
    plt.title('Random Forest Feature Importance (Top 20 Genes)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig('feature_importance_rf.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: feature_importance_rf.png")

    # --- Consensus Top Genes ---
    top50_xgb = set(xgb_imp.head(50)['Gene'])
    top50_rf = set(rf_imp.head(50)['Gene'])
    consensus_genes = sorted(top50_xgb.intersection(top50_rf))

    pd.DataFrame({'Gene': consensus_genes}).to_csv('consensus_top_genes.csv', index=False)
    print("  Saved: consensus_top_genes.csv")

    # ----------------------------------------------------------
    # Final Summary
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED!")
    print("=" * 60)
    print("\nOutput files:")
    print("  - model_comparison_summary.csv")
    print("  - roc_curves_comparison.png")
    print("  - confusion_matrices.png")
    print("  - cv_results_boxplot.png")
    print("  - feature_importance_xgboost.csv/png")
    print("  - feature_importance_rf.csv/png")
    print("  - consensus_top_genes.csv")

    print(f"\nBest Model: {summary_df.iloc[0]['Model']}")
    print(f"  AUC: {summary_df.iloc[0]['AUC']}")

    print(f"\nTop 10 Genes (XGBoost):")
    print(xgb_imp.head(10).to_string(index=False))

    print(f"\nConsensus Top Genes ({len(consensus_genes)} genes in both XGBoost & RF top 50):")
    print(consensus_genes[:20])

    # Best model classification report
    best_name = summary_df.iloc[0]['Model']
    best_res = results[best_name]
    print(f"\n{best_name} Classification Report:")
    print(classification_report(y_test, best_res['y_pred'],
                               target_names=[f'{c}' for c in le.classes_]))


if __name__ == "__main__":
    main()
