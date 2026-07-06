# Ms-dilatometry
Python code to automatically determine Ms from dilatometry curves using the x% transformation method through the lever rule. It is important to know that you need to fully transform to martensite for this method.


The inputs are as follows:

martensite_dilatometry(csv_path,
                        mask_temp=600,
                        FCC_fit_range=(450, 550),
                        BCC_fit_range=(100, 200),
                        smooth=15,
                        target_fraction=0.05,
                        max_Ms_temp=450,
                        Lo = 10000)

csv_path: Full path to the dilatometry .csv file.
mask_temp: Highest temperature displayed in the plot (default 600°C).
FCC_fit_range: Temperature range used for linear fitting of the FCC phase (default 450 to 550 °C).
BCC_fit_range: Temperature range used for linear fitting of the BCC phase (default 100 to 200 °C).
smooth: Smoothing window applied to the derivative curve (default 15 seems to work well).
target_fraction: Fraction transformed used to define Ms (e.g., 0.02 = 2%, 0.05 = 5%, 0.10 = 10%, default is 0.05).
max_Ms_temp: Upper limit for Ms determination (this prevents unrealistically high Ms values that can occasionally appear from poor extrapolation of the FCC fit, particularly when low transformation fractions (e.g., 2%) are used).
Lo: Initial specimen length in μm (default 10000 micrometre).


If you wish to not change any of the variables you can call the function as:
Ms = martensite_dilatometry(r"C:\your_path\...\...\your_file.csv")

or if you wish to change the variables as:

Ms = martensite_dilatometry(r"C:\your_path\...\...\your_file.csv",
                        mask_temp=600,
                        FCC_fit_range=(450, 550),
                        BCC_fit_range=(100, 200),
                        smooth=15,
                        target_fraction=0.05,
                        max_Ms_temp=450,
                        Lo = 10000)

You might have to change the way data is read if you do not use a Bähr DIL 805A/D dilatometer. 
