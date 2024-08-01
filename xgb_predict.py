import joblib
import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import argparse
import os
from pathlib import Path

# python /home/max/stayahead/analysis/scripts/xgb_predict.py -m /home/max/stayahead/analysis/datasets/xgb/outputs/esm/bind_expr/ESM_1step_pred_wuhan_v1_b2b_03-07_model.joblib -x /home/max/stayahead/analysis/datasets/xgb/outputs/esm/bind_expr/X_scaler.joblib -y /home/max/stayahead/analysis/datasets/xgb/outputs/esm/bind_expr/y_scaler.joblib -i /home/max/stayahead/analysis/datasets/current/ESM_BA1_100_b2b_BL.csv -o /home/max/stayahead/analysis/datasets/xgb/outputs/esm_test

def load_scaler(scaler_path):
    return joblib.load(scaler_path)

def preprocess_data(df, feature_columns, x_scaler):
    X = df[feature_columns]
    X_scaled = x_scaler.transform(X)
    return X_scaled

def predict(model, X_scaled, y_scaler):
    dmatrix = xgb.DMatrix(X_scaled)
    y_pred_scaled = model.predict(dmatrix)
    y_pred_rescaled = y_scaler.inverse_transform(y_pred_scaled)
    return y_pred_rescaled

def save_predictions(predictions, output_path):
    # predictions_df = pd.DataFrame(predictions, columns=['bind', 'delta_bind', 'expr', 'delta_expr'])
    predictions_df = pd.DataFrame(predictions, columns=['bind', 'expr'])
    predictions_df.to_csv(output_path, index=False)
    print(f"Predictions saved to {output_path}")

def create_parser():
    parser = argparse.ArgumentParser(
        description='Use XGBoost model to make predictions'
    )
    parser.add_argument(
        "-m", '--model_path', help='Path to the saved model', type=Path, required=True
    )
    parser.add_argument(
        "-x", '--x_scaler_path', help='Path to the saved x_scaler', type=Path, required=True
    )
    parser.add_argument(
        "-y", '--y_scaler_path', help='Path to the saved y_scaler', type=Path, required=True
    )
    parser.add_argument(
        "-i", '--input_file', type=Path, help='Path to the new dataset (CSV file)', required=True
    )
    parser.add_argument(
        "-o", '--output_file', type=Path, help='Path to save the predictions (CSV file)', required=True
    )
    parser.add_argument(
        "-at", "--al_type", help="Algorithm type; alphafold or esmfold", type=str, required=True
    )
    parser.add_argument(
        "-tt", "--target_type", help="Type of target variable; bind_expr, or bind_expr_deltas", type=str, required=True
    )
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()

    model = xgb.XGBRegressor({'n_thread': 4})
    model.load_model(args.model_path)
    x_scaler = load_scaler(args.x_scaler_path)
    y_scaler = load_scaler(args.y_scaler_path)

    df = pd.read_csv(args.input_file)

    # Make sure to use the same feature columns as in the training script
    feature_columns = ['rmsd', 'tm_score', 'sasa', 'avg_hydro', 'plddt', 'local_net_energy','agmata', 'backbone', 'coil', 'disoMine', 'earlyFolding', 'helix', 'ppII', 'sheet', 'sidechain']  # ESM
    X_scaled = preprocess_data(df, feature_columns, x_scaler)

    predictions = predict(model, X_scaled, y_scaler)

    save_predictions(predictions, args.output_file)

if __name__ == "__main__":
    main()