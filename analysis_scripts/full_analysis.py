import joblib
import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, RandomizedSearchCV, GridSearchCV, cross_val_score, KFold
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import argparse
import os
from pathlib import Path

# scikit-learn/1.3.1-gfbf-2023a

# python /home/max/stayahead/analysis/scripts/full_analysis.py -i /home/max/stayahead/analysis/datasets/inputs/5B -o /home/max/stayahead/analysis/outputs/5B
# Define the parameter grid for XGBRegressor
# 5B
# param_dist = {
#     'n_estimators': [800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900],
#     'learning_rate': [0.03, 0.0325, 0.035, 0.0375],
#     'max_depth': [6, 7, 8, 9, 10],
#     'min_child_weight': [5, 6, 7, 8, 9],
#     'gamma': [0, 0.05, 0.075, 0.1],
#     'subsample': [0.6, 0.7, 0.8, 0.9, 1.0],
#     'colsample_bytree': [0.75, 0.8, 0.85],
#     'reg_alpha': [0.15, 0.175, 0.2, 0.225, 0.25],
#     'reg_lambda': [1.5, 1.75, 2, 2.25, 2.5]
# }

# 16B
# param_dist = {
#     'n_estimators': [800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900],
#     'learning_rate': [0.0275, 0.03, 0.0325, 0.035, 0.0375],
#     'max_depth': [5, 6, 7, 8, 9],
#     'min_child_weight': [4, 5, 6, 7, 8],
#     'gamma': [0, 0.025, 0.05, 0.075],
#     'subsample': [0.6, 0.7, 0.8, 0.9, 1.0],
#     'colsample_bytree': [0.75, 0.8, 0.85],
#     'reg_alpha': [0, 0.025, 0.05, 0.075, 0.1],
#     'reg_lambda': [0.5, 0.75, 1, 1.25, 1.5]
# }

# 1H-6B
# param_dist = {
#     'n_estimators': [800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900],
#     'learning_rate': [0.025, 0.0275, 0.03, 0.0325, 0.035],
#     'max_depth': [7, 8, 9, 10, 11],
#     'min_child_weight': [4, 5, 6, 7, 8],
#     'gamma': [0, 0.01, 0.02, 0.03],
#     'subsample': [0.6, 0.7, 0.8, 0.9, 1.0],
#     'colsample_bytree': [0.75, 0.8, 0.85],
#     'reg_alpha': [0.05, 0.075, 0.1, 0.15, 0.175],
#     'reg_lambda': [1.5, 1.75, 2, 2.25, 2.5]
# }

# 1H-20B
param_dist = {
    'n_estimators': [900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000],
    'learning_rate': [0.03, 0.0325, 0.035, 0.0375, 0.04],
    'max_depth': [6, 7, 8, 9, 10],
    'min_child_weight': [6, 7, 8, 9, 10],
    'gamma': [0.025, 0.05, 0.075, 0.1, 0.125],
    'subsample': [0.6, 0.7, 0.8, 0.9, 1.0],
    'colsample_bytree': [0.75, 0.8, 0.85],
    'reg_alpha': [0, 0.025, 0.05, 0.075, 0.1],
    'reg_lambda': [0.5, 0.75, 1, 1.25, 1.5]
}

def gboost_regressor(df, param_dist, random_search, model_save_path, log_file_path, plot_save_path):
    # Preprocessing
    df = df.dropna(subset=['bind', 'delta_bind', 'expr', 'delta_expr', 'confidence_bind', 'confidence_expr'])
    df['pi_score'] = pd.to_numeric(df['pi_score'], errors='coerce')

    # Separate features and target variables
    # X = df.drop(columns=['seq_id', 'wildtype', 'site', 'mutation', 'bind', 'delta_bind', 'expr', 'delta_expr', 'confidence_bind', 'confidence_expr', 'sb', 'global_net_energy', 'total_hydro', 'Num_intf_residues', 'contact_pairs', 'hb', 'Hydrophobhic', 'Polar', 'Charged', 'iptm', 'iptm_ptm'])
    X = df.drop(columns=['seq_id', 'wildtype', 'mutation', 'bind', 'delta_bind', 'expr', 'delta_expr', 'confidence_bind', 'confidence_expr', 'global_net_energy', 'total_hydro'])
    y = df[['bind', 'delta_bind', 'expr', 'delta_expr']]

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Standardization
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Scale the target variables
    y_scaler = StandardScaler()
    y_train_scaled = y_scaler.fit_transform(y_train)
    y_test_scaled = y_scaler.transform(y_test)

    # Initialize XGBoost model
    model = xgb.XGBRegressor(objective='reg:squarederror', random_state=42)

    log_lines = []

    if random_search:
        # Initialize RandomizedSearchCV for the base XGBRegressor
        random_search = RandomizedSearchCV(estimator=model, param_distributions=param_dist, 
                                        n_iter=600, scoring='neg_mean_squared_error', 
                                        cv=5, verbose=2, random_state=42, n_jobs=-1)

        # Fit the random search model for the base XGBRegressor
        random_search.fit(X_train_scaled, y_train_scaled[:, 0])  # Fit using the first target variable to get the best params

        # Log the best parameters and best score
        log_lines.append(f'Best Parameters: {random_search.best_params_}')
        log_lines.append(f'Best Score: {random_search.best_score_}')

        # Train the MultiOutputRegressor with the best parameters
        best_params = random_search.best_params_
    else:
        # Initialize GridSearchCV
        grid_search = GridSearchCV(estimator=model, param_grid=param_dist, 
                                   scoring='neg_mean_squared_error', 
                                   cv=5, verbose=2, n_jobs=-1)

        # Fit the grid search model
        grid_search.fit(X_train_scaled, y_train_scaled[:, 0])

        # Log the best parameters and best score
        log_lines.append(f'Best Parameters: {grid_search.best_params_}')
        log_lines.append(f'Best Score: {grid_search.best_score_}')

        # Train the model with best parameters
        best_params = grid_search.best_params_
        
    base_model = xgb.XGBRegressor(**best_params, objective='reg:squarederror', random_state=42)
    multi_output_model = MultiOutputRegressor(base_model)

    # Cross-validation setup
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_results = cross_val_score(multi_output_model, X_train_scaled, y_train_scaled, cv=kf, scoring='neg_mean_squared_error')

    # Fit the multi-output model
    multi_output_model.fit(X_train_scaled, y_train_scaled)

    # Save the model
    joblib.dump(multi_output_model, model_save_path)
    log_lines.append(f'Model saved to {model_save_path}')

    # Make predictions
    y_pred_scaled = multi_output_model.predict(X_test_scaled)

    # Rescale predictions back to original scale
    y_pred_rescaled = y_scaler.inverse_transform(y_pred_scaled)

    # Evaluate the model for each target variable
    for i, col in enumerate(y.columns):
        mse = mean_squared_error(y_test.iloc[:, i], y_pred_rescaled[:, i])
        r2 = r2_score(y_test.iloc[:, i], y_pred_rescaled[:, i])
        log_lines.append(f'{col} - Mean Squared Error: {mse}')
        log_lines.append(f'{col} - R-squared: {r2}')

    # Log cross-validation results
    log_lines.append(f'Cross-Validation MSE Scores: {-cv_results}')
    log_lines.append(f'Cross-Validation Mean MSE: {-cv_results.mean()}')
    log_lines.append(f'Cross-Validation Standard Deviation of MSE: {cv_results.std()}')

    # Write log lines to a file
    with open(log_file_path, 'w') as log_file:
        for line in log_lines:
            log_file.write(line + '\n')

    # Feature importance (using the first model in the MultiOutputRegressor)
    importance = multi_output_model.estimators_[0].get_booster().get_score(importance_type='weight')
    importance = {X.columns[int(k[1:])]: v for k, v in importance.items()}
    sorted_importance = dict(sorted(importance.items(), key=lambda item: item[1], reverse=True))

    plt.figure(figsize=(10, 8))
    plt.barh(range(len(sorted_importance)), list(sorted_importance.values()), align='center')
    plt.yticks(range(len(sorted_importance)), list(sorted_importance.keys()))
    plt.xlabel('F score')
    plt.ylabel('Features')
    plt.title('Feature importance')
    plt.savefig(plot_save_path)  # Save the plot to a file
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
        "-r", '--random', type=str, help='Activate random GridSearchCV', required=False
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
        df = pd.read_csv(args.input_dir / file)
        name = file.split('.')[0]
        model_save_path = os.path.join(args.output_path, f'{name}_model.joblib')
        plot_save_path = os.path.join(args.output_path, f'{name}_feature_importance.png')
        log_file_path = os.path.join(args.output_path, f'{name}_log.txt')
        gboost_regressor(df, param_dist, args.random, model_save_path, log_file_path, plot_save_path)

if __name__ == "__main__":
    main()