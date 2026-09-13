# Electric Field and Equipotential Visualizer

A computational electrostatics tool for modeling and analyzing the electric field and electric potential produced by multiple point charges.

The simulator applies **Coulomb's law and the principle of superposition** across a configurable 2D spatial grid to generate field streamlines, equipotential contours, and field-magnitude maps.

## Key Features

* Custom and predefined point-charge configurations
* Electric field and equipotential visualization
* Electric-field magnitude mapping
* Numerical validation of

  $$
  \mathbf{E} \approx -\nabla V
  $$
* Grid-resolution convergence analysis
* Probe-point field and potential calculations
* CSV export of computed field data

## Built With

**Python · NumPy · Matplotlib**

The project demonstrates how fundamental electrostatic laws can be translated into a computational model and **quantitatively validated through numerical analysis**.

## Run

```bash
pip install numpy matplotlib
python electric_field_visualizer.py
```

Generated figures and numerical data are saved to the `outputs/` directory.

**Author:** Durodola Daniel

## Simulations
<img width="1282" height="603" alt="field_validation" src="https://github.com/user-attachments/assets/3940cf65-fe4d-4f92-8792-575c886dee14" />

<img width="1282" height="614" alt="field_orthogonality" src="https://github.com/user-attachments/assets/88d429ac-c05b-4365-8966-9465b8ef0083" />

<img width="1282" height="577" alt="field_equipotential" src="https://github.com/user-attachments/assets/0ee3228d-24a5-4d35-836f-da613d254f77" />

