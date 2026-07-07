import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
import os

def martensite_dilatometry(csv_path,
                        mask_temp=600,
                        FCC_fit_range=(450, 550),
                        BCC_fit_range=(100, 200),
                        smooth=15,
                        target_fraction=0.05,
                        max_Ms_temp=450,
                        Lo = 10000):


    # --- Load CSV header ---
    with open(csv_path, encoding='unicode_escape') as fh:
        reader = csv.reader(fh, delimiter=';')
        next(reader)
        next(reader)
        next(reader)
        next(reader)
        csv_header = next(reader)

    # --- Load full dataset ---
    df = pd.read_csv(csv_path, delimiter=';', encoding='unicode_escape',
                     header=None, skiprows=5, names=csv_header)

    # --- Cooling detection ---
    T_raw = df['Temperature'].values
    T_smooth = pd.Series(T_raw).rolling(50, center=True, min_periods=1).mean()
    dTdt = np.gradient(T_smooth)

    N = 20
    cool_start = None
    for i in range(len(dTdt) - N):
        if np.all(dTdt[i:i+N] < 0):
            cool_start = i
            break

    print("Cooling starts at index:", cool_start)

    # --- Apply cooling cutoff ---
    T_cool = T_raw[cool_start:]
    dL_cool = df['Change in Length'].values[cool_start:]
    t_cool = df['Time'].values[cool_start:]

    # --- Mask below chosen temperature ---
    mask = T_cool < mask_temp
    T = T_cool[mask]
    dL = dL_cool[mask]
    t = t_cool[mask]

    # --- Normalize ---
    dL_fit = dL * (100 / Lo)

    # --- Fit FCC/BCC lines ---
    FCC_mask = (T >= FCC_fit_range[0]) & (T <= FCC_fit_range[1])
    BCC_mask = (T >= BCC_fit_range[0]) & (T <= BCC_fit_range[1])

    FCC_coeff = np.polyfit(T[FCC_mask], dL_fit[FCC_mask], 1)
    BCC_coeff = np.polyfit(T[BCC_mask], dL_fit[BCC_mask], 1)

    FCC_fit = np.poly1d(FCC_coeff)
    BCC_fit = np.poly1d(BCC_coeff)

    FCC_curve = FCC_fit(T)
    BCC_curve = BCC_fit(T)

    # --- Martensite fraction ---
    BCC_fraction = (dL_fit - FCC_curve) / (BCC_curve - FCC_curve)

    # --- FIND TEMPERATURE WHERE BCC FRACTION = x ---
    target = target_fraction
    max_valid_T = max_Ms_temp

    # Use the SAME arrays used for plotting
    T_arr = np.asarray(T)
    f_arr = np.asarray(BCC_fraction)

    # Difference from target
    diff = f_arr - target

    # Sign-change crossing detection
    cross_idx = np.where(diff[:-1] * diff[1:] <= 0)[0]

    # Only accept crossings below max_valid_T
    cross_idx = [i for i in cross_idx if T_arr[i] < max_valid_T]

    if not cross_idx:
        print(f"No intersection with f={target} below {max_valid_T} °C")
        T_intersect = None
    else:
        i = cross_idx[0]

        # Linear interpolation between the two points that bracket the crossing
        x0, x1 = T_arr[i], T_arr[i+1]
        y0, y1 = f_arr[i], f_arr[i+1]

        T_intersect = x0 + (target - y0) * (x1 - x0) / (y1 - y0)
        print(f"Intersection at f={target} occurs at T = {T_intersect:.2f} °C")


    # --- Derivative ---
    derivative_curve = []
    for i in range(smooth, len(BCC_fraction)):
        dfdT = -(BCC_fraction[i] - BCC_fraction[i-smooth]) / (T[i] - T[i-smooth])
        derivative_curve.append(dfdT)

    # --- Plotting ---
    fig, ax = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Length change
    ax[0,0].set_title('Length change during cooling')
    ax[0,0].plot(T, dL_fit, linewidth=3)
    ax[0,0].set_xlabel('Temperature [°C]')
    ax[0,0].set_ylabel('ΔL/L₀ [%]')
    ax[0,0].grid()

    # 2. FCC/BCC fits
    ax[0,1].set_title('FCC/BCC linear fits')
    ax[0,1].plot(T, dL_fit, label='Measured')
    ax[0,1].plot(T, FCC_curve, label='FCC fit')
    ax[0,1].plot(T, BCC_curve, label='BCC fit')
    ax[0,1].set_xlabel('Temperature [°C]')
    ax[0,1].set_ylabel('ΔL/L₀ [%]')
    ax[0,1].grid()
    ax[0,1].legend()

    # 3. BCC fraction
    ax[1,0].set_title('Martensite phase fraction')
    ax[1,0].plot(T, BCC_fraction, linewidth=3)
    ax[1,0].axhline(target_fraction, color='red', linestyle='--', linewidth=2)
    if T_intersect is not None:
        ax[1,0].plot(T_intersect, target_fraction, 'ro')
        ax[1,0].annotate(f"{T_intersect:.0f} °C",
                         (T_intersect, target_fraction),
                         textcoords="offset points", xytext=(10,10))
    ax[1,0].set_xlabel('Temperature [°C]')
    ax[1,0].set_ylabel('Martensite fraction')
    ax[1,0].grid()

    # 4. Derivative
    ax[1,1].set_title('df/dT')
    ax[1,1].plot(T[smooth:], derivative_curve)
    ax[1,1].set_xlabel('Temperature [°C]')
    ax[1,1].set_ylabel('ΔBCC/ΔT')
    ax[1,1].grid()

    sample_name = os.path.basename(csv_path).replace(".csv", "")
    fig.suptitle(f"Dilatometry Analysis – {sample_name}", fontsize=15, y=1.005)
    plt.tight_layout()
    plt.show()

    return T_intersect
