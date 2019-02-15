import numpy as np
import matplotlib.pyplot as plt
import load_inputs as curves

if __name__ == "__main__":
    # mode1_composite

    # left panel
    mode1_composite_fig, mode1_composite_axes = plt.subplots(nrows=1, ncols=2)
    mode1_composite_fig.set_size_inches(10, 5)
    mode1_composite_mean_N_Black_x = np.transpose(curves.mode1_composite_mean_N_Black[1:])[0]
    mode1_composite_mean_N_Black_y = np.transpose(curves.mode1_composite_mean_N_Black[1:])[1]
    mode1_composite_Err_N_Black_Upper_x = np.transpose(curves.mode1_composite_Err_N_Black_Upper[1:])[0]
    mode1_composite_Err_N_Black_Upper_y = np.transpose(curves.mode1_composite_Err_N_Black_Upper[1:])[1]
    mode1_composite_Err_N_Black_Lower_x = np.transpose(curves.mode1_composite_Err_N_Black_Lower[1:])[0]
    mode1_composite_Err_N_Black_Lower_y = np.transpose(curves.mode1_composite_Err_N_Black_Lower[1:])[1]

    mode1_composite_mean_N_Red_x = np.transpose(curves.mode1_composite_mean_N_Red[1:])[0]
    mode1_composite_mean_N_Red_y = np.transpose(curves.mode1_composite_mean_N_Red[1:])[1]
    mode1_composite_Err_N_Red_Upper_x = np.transpose(curves.mode1_composite_Err_N_Red_Upper[1:])[0]
    mode1_composite_Err_N_Red_Upper_y = np.transpose(curves.mode1_composite_Err_N_Red_Upper[1:])[1]
    mode1_composite_Err_N_Red_Lower_x = np.transpose(curves.mode1_composite_Err_N_Red_Lower[1:])[0]
    mode1_composite_Err_N_Red_Lower_y = np.transpose(curves.mode1_composite_Err_N_Red_Lower[1:])[1]

    mode1_composite_fig_x_ticks = ["%.2f" % el for el in np.transpose(curves.mode1_composite_mean_N_Black[1::int(len(curves.mode1_composite_mean_N_Black)/8)])[0]]
    mode1_composite_fig_y_ticks = [0.2 * r for r in range(6)]

    mode1_composite_axes[0].plot(mode1_composite_mean_N_Black_x, mode1_composite_mean_N_Black_y, 'k')
    mode1_composite_axes[0].plot(mode1_composite_Err_N_Black_Upper_x, mode1_composite_Err_N_Black_Upper_y, 'k--')
    mode1_composite_axes[0].plot(mode1_composite_Err_N_Black_Lower_x, mode1_composite_Err_N_Black_Lower_y, 'k--')
    mode1_composite_axes[0].plot(mode1_composite_mean_N_Red_x, mode1_composite_mean_N_Red_y, 'r')
    mode1_composite_axes[0].plot(mode1_composite_Err_N_Red_Upper_x, mode1_composite_Err_N_Red_Upper_y, 'r--')
    mode1_composite_axes[0].plot(mode1_composite_Err_N_Red_Lower_x, mode1_composite_Err_N_Red_Lower_y, 'r--')
    mode1_composite_axes[0].set_xlabel('c')
    mode1_composite_axes[0].set_ylabel('<n>/kpt')

    # right panel
    mode1_error_compare_heuristic_x = np.transpose(curves.mode1_error_compare_heuristic[1:])[0]
    mode1_error_compare_heuristic_y = np.transpose(curves.mode1_error_compare_heuristic[1:])[1]

    mode1_composite_axes[1].plot(mode1_error_compare_heuristic_x, mode1_error_compare_heuristic_y, 'k')
    mode1_composite_axes[1].set_xlabel('x')
    mode1_composite_axes[1].set_ylabel(r'$\delta x^{2}$/$x^{2}$')
    mode1_composite_axes[1].set_ylim([0, 1.5])

    plt.show()

    # mode1_MLE
    #mode1_composite_fig_x_values = np.transpose(curves.mode1_MLE_compare_heuristic[1:])[0]
    #mode1_composite_fig_y_values = np.transpose(curves.mode1_MLE_compare_heuristic[1:])[1]
    #mode1_composite_fig_x_ticks = np.transpose(curves.mode1_MLE_compare_heuristic[1::40])[0]
    #mode1_composite_fig_y_ticks = np.transpose(curves.mode1_MLE_compare_heuristic[1::40])[1]
    #plt.xticks(mode1_composite_fig_x_ticks, ["%.2f" % el for el in mode1_composite_fig_x_ticks])
    #mode1_composite_axes[0].plot(mode1_composite_fig_x_values[0:399], mode1_composite_fig_y_values[0:399], 'r')
    #mode1_composite_axes[0].plot(mode1_composite_fig_x_values[400:], mode1_composite_fig_y_values[400:], 'r')
    #plt.show()

