import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

def plot_sequins_ladder(
    df, x_col, y_col, pseudocount=0,
    x_label=None, y_label=None, 
    title=None, filename=None
):
    """
    Plot a log-log regression of a samples Sequins ladder based on known input concentration vs normalised abundance.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing input concentrations and normalised observed abundance.
    x_col : str
        Column name for input concentrations (e.g., gcul, ngul).
    y_col : str
        Column name for measured abundance (e.g., TPM, CPM).
    pseudocount : float, default 0.0
        Pseudocount added to y values before log transformation to avoid log(0).
    x_label : str, optional
        Label for the x-axis. If None, auto-generates from x_col.
    y_label : str, optional
        Label for the y-axis. If None, auto-generates from y_col.
    title : str, optional
        Title for the plot. If None, auto-generates from axis labels.
    filename : str, optional
        File path to save the plot. If None, generates a filename from column names.

    Returns
    -------
    statsmodels.regression.linear_model.RegressionResultsWrapper
        Fitted OLS regression model.
    """
    required_cols = {x_col, y_col}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")

    df = df.copy()
    df = df.sort_values(by=x_col)

    # Fit regression model on log10-transformed data
    X = sm.add_constant(np.log10(df[x_col]))
    
    y_log = np.log10(df[y_col]+ pseudocount)
    model = sm.OLS(y_log, X).fit()
    df['y_pred_log'] = model.predict(X)

    # Extract regression statistics
    slope = model.params.iloc[1]
    r_squared = model.rsquared

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(df[x_col], df[y_col]+pseudocount,
               marker='o', facecolors='none', edgecolors='blue', s=70, label='Sequins points')
    ax.plot(df[x_col], 10 ** df['y_pred_log'],
            color='black', label='Sequins Ladder Regression fit', linewidth=1)

    ax.set_xscale('log')
    ax.set_yscale('log')

    # Axis labels and title
    x_label = x_label or f'Input ({x_col})'
    y_label = y_label or f'Abundance ({y_col})'
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)

    title = title or f'Sequins Ladder - {y_label} vs {x_label}'
    ax.set_title(title, fontsize=14)

     # R² and slope annotation on plot
    ax.text(0.95, 0.05,
            f"$R^2$: {r_squared:.4f}\nSlope: {slope:.4f}",
            fontsize=12, color='black', transform=ax.transAxes, ha='right',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85))

    ax.legend(fontsize=12, loc='upper left')
    ax.grid(True, which='major', linestyle='--', linewidth=0.5)

    # File saving
    filename = filename or f"Sequins_{y_col}_vs_{x_col}_loglog.png"
    os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.show()

    return model


def limits_calculation(
    model: sm.regression.linear_model.RegressionResultsWrapper,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    cv_threshold: float = 35.0,
    pseudocount: float = 0.0
) -> (float, float, pd.DataFrame):
    """
    Calculate the Limit of Quantification (LoQ) and Limit of Detection (LoD)
    from a log–log OLS model of the Sequins metagenomics ladder.

    Parameters
    ----------
    model : statsmodels RegressionResultsWrapper
        A fitted OLS model from log10(y) ~ const + log10(x).
    df : pd.DataFrame
        The original DataFrame that was used to fit the model.
    x_col : str
        Name of the column in `df` containing the known input concentration (x) values.
    y_col : str
        Name of the column in `df` containing the observed normalised abudannce (y) values.
    cv_threshold : float, default 35.0
        Maximum CV% allowed for back calculated concentration values at a given known input concentration (x) value (for LoQ).
    pseudocount : float, default 0.0
        Pseudocount added to y values before log transformation to avoid log(0).

    Returns
    -------
    loq : float
        The lowest known input concentration (x) value at which all sequins were >0
        and CV% of back-calculated input concentrations <= cv_threshold.
    lod : float
        The lowest known input concentration (x) at which all sequins were >0.
    results_df : pd.DataFrame
        One row per unique known input concentration (x) with columns:
          - known_input_concentration
          - all_detected
          - mean_backcalc
          - std_backcalc
          - cv_percent
    """
    if x_col not in df.columns or y_col not in df.columns:
        raise ValueError(f"Columns '{x_col}' and/or '{y_col}' not found in the DataFrame.")

    # extract intercept and slope from the model
    b0, b1 = model.params

    # compute back‑calculated input concentrations (x values) for each abundance (y) value
    tmp = df[[x_col, y_col]].copy()
    tmp['log_y'] = np.log10(tmp[y_col]+ pseudocount)
    tmp['log_x_back'] = (tmp['log_y'] - b0) / b1
    tmp['x_back'] = 10 ** tmp['log_x_back']

    # summarise by known input concentrations (x values)
    records = []
    for known, group in tmp.groupby(x_col):
        obs_vals = group[y_col].values
        all_detected = np.all(obs_vals > 0)
        back = group['x_back'].values
        mean_bc = np.mean(back)
        std_bc  = np.std(back, ddof=1)
        cv = (std_bc / mean_bc * 100) if mean_bc > 0 else np.nan

        records.append({
            'known_input_concentration':     known,
            'all_detected':  all_detected,
            'mean_backcalc': mean_bc,
            'std_backcalc':  std_bc,
            'cv_percent':    cv
        })

    results_df = pd.DataFrame(records).sort_values(by='known_input_concentration')
    # Compute LoD: smallest nominal_x where all_detected is True
    lod_candidates = results_df[results_df['all_detected']]
    lod = lod_candidates['known_input_concentration'].min() if not lod_candidates.empty else np.nan

    # Compute LoQ: smallest nominal_x where all_detected & cv ≤ threshold
    loq_candidates = lod_candidates[lod_candidates['cv_percent'] <= cv_threshold]
    loq = loq_candidates['known_input_concentration'].min() if not loq_candidates.empty else np.nan

    results_df['known_input_concentration'] = results_df['known_input_concentration'].apply(lambda x: f'{x:.3e}')
    results_df['mean_backcalc'] = results_df['mean_backcalc'].apply(lambda x: f'{x:.3e}')
    results_df['std_backcalc'] = results_df['std_backcalc'].apply(lambda x: f'{x:.3e}')

    return loq, lod, results_df