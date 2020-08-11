# cellsignalling
Code for the publication *Effects of cross-talk and pleiotropy on the specificity and accuracy of receptor signaling*

All figures can be generated from the script */figures/generate_paper_figures.py*

- Figure 1D: run the function `figure_1_and_2_heatmaps()`; saves to /output/ligand1
- Figure 1E&2B: run the function `plot_1E_and_2B()`; saves to /output
- Figure 1F: run the function `plot_1F()`; saves to /output
- Figure 2C&2D: run the function `figure_1_and_2_heatmaps()`; saves to /output/ligand1
- Figure 3C&D: run the function `multiple_heatmaps()`; saves to /output/ligand2
- Figure S1: run the function `plot_S1()`; saves to /output
- Figure S2: run the function `plot_S2()`; saves to /output
- Figure S3: run the function `plot_S3()`; saves to /output/ligand1/SI_ratios
- Figure S4: run *make_supplementary_histograms.py* in /simulation; saves to /simulation/output

To recreate all figures:
1. Run the Mathematica master file in the parent directory
2. run */figures/equations_txt2python.py*
3. run */figures/generate_paper_figures.py*
4. run */simulation/make_supplementary_histograms.py*
