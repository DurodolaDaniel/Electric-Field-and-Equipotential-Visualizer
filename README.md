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

## Example Results

![Electric Field and Equipotential](outputs/field_equipotential.png)

![Electric Field Magnitude](outputs/electric_field_magnitude.png)

![Numerical Validation](outputs/field_validation.png)

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
