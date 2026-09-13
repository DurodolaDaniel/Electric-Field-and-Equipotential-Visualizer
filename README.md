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
<img width="991" height="603" alt="WhatsApp Image 2026-09-13 at 14 20 31" src="https://github.com/user-attachments/assets/f2e3bca3-c084-435c-bb23-dcefd99bc4d0" />


<img width="991" height="614" alt="WhatsApp Image 2026-09-13 at 14 20 31 (1)" src="https://github.com/user-attachments/assets/dedbda79-08d9-4366-b5fe-a22a424ef418" />


<img width="1282" height="577" alt="field_equipotential" src="https://github.com/user-attachments/assets/0564567e-b507-4644-ae0f-27548fbf47de" />

