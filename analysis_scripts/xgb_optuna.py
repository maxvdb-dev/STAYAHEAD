import joblib
import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, KFold
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import argparse
import os
from pathlib import Path
import optuna
import json
# import warnings

#python /Users/maxvandenboom/stayahead/analysis/scripts/xgb_optuna.py -i /Users/maxvandenboom/stayahead/analysis/datasets/inputs/1H_6B -o /Users/maxvandenboom/stayahead/analysis/outputs/1H_6B -p /Users/maxvandenboom/stayahead/analysis/datasets/inputs/1H_6B/params.json 

# Suppress specific warnings
# warnings.filterwarnings('ignore', category=FutureWarning, module='sklearn')

def load_param_ranges(json_path):
    with open(json_path, 'r') as f:
        param_ranges = json.load(f)
    return param_ranges

def objective(trial, X_train_scaled, y_train_scaled, param_ranges):
    param = {
        'verbosity': 0,
        'objective': 'reg:squarederror',
        'booster': trial.suggest_categorical('booster', param_ranges['booster']),
        'lambda': trial.suggest_float('lambda', *param_ranges['lambda']),
        'alpha': trial.suggest_float('alpha', *param_ranges['alpha']),
        'subsample': trial.suggest_float('subsample', *param_ranges['subsample']),
        'colsample_bytree': trial.suggest_float('colsample_bytree', *param_ranges['colsample_bytree']),
    }

    if param['booster'] in ['gbtree', 'dart']:
        param['max_depth'] = trial.suggest_int('max_depth', *param_ranges['max_depth'])
        param['min_child_weight'] = trial.suggest_int('min_child_weight', *param_ranges['min_child_weight'])
        param['eta'] = trial.suggest_float('eta', *param_ranges['eta'])
        param['gamma'] = trial.suggest_float('gamma', *param_ranges['gamma'])
        param['grow_policy'] = trial.suggest_categorical('grow_policy', param_ranges['grow_policy'])

    if param['booster'] == 'dart':
        param['sample_type'] = trial.suggest_categorical('sample_type', param_ranges['sample_type'])
        param['normalize_type'] = trial.suggest_categorical('normalize_type', param_ranges['normalize_type'])
        param['rate_drop'] = trial.suggest_float('rate_drop', *param_ranges['rate_drop'])
        param['skip_drop'] = trial.suggest_float('skip_drop', *param_ranges['skip_drop'])

    model = xgb.XGBRegressor(**param)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    rmse_scores = []
    for train_index, valid_index in kf.split(X_train_scaled):
        X_tr, X_val = X_train_scaled[train_index], X_train_scaled[valid_index]
        y_tr, y_val = y_train_scaled[train_index], y_train_scaled[valid_index]

        dtrain = xgb.DMatrix(X_tr, label=y_tr)
        dvalid = xgb.DMatrix(X_val, label=y_val)

        model.fit(X_tr, y_tr)
        y_pred = model.predict(X_val)
        rmse = mean_squared_error(y_val, y_pred, squared=False)
        rmse_scores.append(rmse)

    return np.mean(rmse_scores)

def gboost_regressor(df, param_ranges, model_save_path, log_file_path, plot_save_path, scaler_save_path, target_type, al_type, dump_path):
    # Preprocessing
    df = df.dropna(subset=['bind', 'delta_bind', 'expr', 'delta_expr', 'confidence_bind', 'confidence_expr'])

    # X = df.drop(columns=['seq_id', 'wildtype', 'site', 'mutation', 'bind', 'delta_bind', 'expr', 'delta_expr', 'confidence_bind', 'confidence_expr', 'global_net_energy', 'total_hydro'])
    X = df.drop(columns=['seq_id', 'wildtype', 'site', 'mutation', 'bind', 'delta_bind', 'expr', 'delta_expr', 'confidence_bind', 'confidence_expr', 'global_net_energy', 'total_hydro', 'sb', 'contact_pairs', 'Num_intf_residues', 'iptm', 'hb', 'iptm_ptm', 'Charged', 'sc', 'Hydrophobhic'])
    if target_type == 'bind_expr_deltas':
        y = df[['bind', 'delta_bind', 'expr', 'delta_expr']]
    elif target_type == 'bind_expr':
        y = df[['bind', 'expr']]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    y_scaler = StandardScaler()
    y_train_scaled = y_scaler.fit_transform(y_train)
    y_test_scaled = y_scaler.transform(y_test)

    study = optuna.create_study(direction='minimize')
    study.optimize(lambda trial: objective(trial, X_train_scaled, y_train_scaled, param_ranges), n_trials=2500, timeout=600)

    best_params = study.best_params
    log_lines = []
    log_lines.append(f'Best Parameters: {best_params}')
    
    base_model = xgb.XGBRegressor(**best_params, objective='reg:squarederror', random_state=42)
    base_model.fit(X_train_scaled, y_train_scaled)

    # joblib.dump(base_model, model_save_path)
    # log_lines.append(f'Model saved to {model_save_path}')

    base_model.save_model(model_save_path)
    # base_model.dump_model(f'{dump_path}_dump.raw.txt', f'{dump_path}_featmap.txt')

    # Save the scalers
    joblib.dump(scaler, os.path.join(scaler_save_path, 'X_scaler.joblib'))
    joblib.dump(y_scaler, os.path.join(scaler_save_path, 'y_scaler.joblib'))
    log_lines.append(f'Scalers saved to {scaler_save_path}')

    y_pred_scaled = base_model.predict(X_test_scaled)
    y_pred_rescaled = y_scaler.inverse_transform(y_pred_scaled)

    for i, col in enumerate(y.columns):
        mse = mean_squared_error(y_test.iloc[:, i], y_pred_rescaled[:, i])
        r2 = r2_score(y_test.iloc[:, i], y_pred_rescaled[:, i])
        log_lines.append(f'{col} - Mean Squared Error: {mse}')
        log_lines.append(f'{col} - R-squared: {r2}')

    with open(log_file_path, 'w') as log_file:
        for line in log_lines:
            log_file.write(line + '\n')

    importance = base_model.get_booster().get_score(importance_type='weight')
    importance = {X.columns[int(k[1:])]: v for k, v in importance.items()}
    sorted_importance = dict(sorted(importance.items(), key=lambda item: item[1], reverse=True))

    plt.figure(figsize=(10, 8))
    plt.barh(range(len(sorted_importance)), list(sorted_importance.values()), align='center')
    plt.yticks(range(len(sorted_importance)), list(sorted_importance.keys()))
    plt.xlabel('F score')
    plt.ylabel('Features')
    plt.title('Feature importance')
    plt.savefig(plot_save_path)
    plt.close()

def create_parser():
    parser = argparse.ArgumentParser(
        description='Train XGBoost model'
    )
    parser.add_argument(
        "-i", '--input_dir', help='Path to the data directory', type=Path, required=True
    )
    parser.add_argument(
        "-o", '--output_path', type=str, help='Path to the output directory', required=True
        )
    parser.add_argument(
        "-p", '--param_json', type=str, help='Path to the hyperparameter JSON file', required=True
    )
    parser.add_argument(
        "-tt", "--target_type", help="Type of target variable; bind_expr, or bind_expr_deltas", type=str, required=True
    )
    parser.add_argument(
        "-at", "--al_type", help="Algorithm type; alphafold or esmfold", type=str, required=True
    )
    return parser

def list_files(directory):
    try:
        files = os.listdir(directory)
        files = [f for f in files if os.path.isfile(os.path.join(directory, f))]
        return files
    except FileNotFoundError:
        print(f"The directory {directory} does not exist.")
        return []

def main():
    parser = create_parser()
    args = parser.parse_args()

    files = list_files(args.input_dir)
    for file in files:
        if file.endswith('.csv'):
            df = pd.read_csv(args.input_dir / file)
        else:
            continue
        name = file.split('.')[0]
        model_save_path = os.path.join(args.output_path, f'{name}_model.json')
        plot_save_path = os.path.join(args.output_path, f'{name}_feature_importance.png')
        log_file_path = os.path.join(args.output_path, f'{name}_log.txt')
        scaler_save_path = args.output_path
        dump_path = args.output_path
        target_type = args.target_type
        al_type = args.al_type
        param_ranges = load_param_ranges(args.param_json)
        gboost_regressor(df, param_ranges, model_save_path, log_file_path, plot_save_path, scaler_save_path, target_type, al_type, dump_path)

if __name__ == "__main__":
    main()
