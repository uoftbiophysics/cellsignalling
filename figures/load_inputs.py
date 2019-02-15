import csv
import os

def load_input(filename):
    with open(os.path.join(os.getcwd(), 'input', filename + '.csv'), 'rb') as f:
        reader = csv.reader(f)
        read_list = list(reader)
        numeric_read_list = [(read_list[0][0],read_list[0][1])] + [(float(el[0]), float(el[1])) for el in read_list[1:]]
    return numeric_read_list

combined_error_composite_compare_c_fullFisher = load_input('combined_error_composite_compare_c_fullFisher')
combined_error_composite_compare_c_heuristic = load_input('combined_error_composite_compare_c_heuristic')
combined_error_composite_compare_c_saddlePointFisher = load_input('combined_error_composite_compare_c_saddlePointFisher')
combined_error_composite_compare_koff_fullFisher = load_input('combined_error_composite_compare_koff_fullFisher')
combined_error_composite_compare_koff_heuristic = load_input('combined_error_composite_compare_koff_heuristic')
combined_error_composite_compare_koff_saddlePointFisher = load_input('combined_error_composite_compare_koff_saddlePointFisher')
combined_MLE_composite_compare_c_heuristic = load_input('combined_MLE_composite_compare_c_heuristic')
combined_MLE_composite_compare_c_numeric = load_input('combined_MLE_composite_compare_c_numeric')
combined_MLE_composite_compare_koff_heuristic = load_input('combined_MLE_composite_compare_koff_heuristic')
combined_MLE_composite_compare_koff_numeric = load_input('combined_MLE_composite_compare_koff_numeric')

KPR_error_composite_compare_c_heuristic = load_input('KPR_error_composite_compare_c_heuristic')
KPR_error_composite_compare_c_saddlePointFisher = load_input('KPR_error_composite_compare_c_saddlePointFisher')
KPR_error_composite_compare_koff_heuristic = load_input('KPR_error_composite_compare_koff_heuristic')
KPR_error_composite_compare_koff_saddlePointFisher = load_input('KPR_error_composite_compare_koff_saddlePointFisher')
KPR_MLE_composite_compare_c_heuristic = load_input('KPR_MLE_composite_compare_c_heuristic')
KPR_MLE_composite_compare_c_numeric = load_input('KPR_MLE_composite_compare_c_numeric')
KPR_MLE_composite_compare_koff_heuristic = load_input('KPR_MLE_composite_compare_koff_heuristic')
KPR_MLE_composite_compare_koff_numeric = load_input('KPR_MLE_composite_compare_koff_numeric')

mode1_error_compare_fullFisher = load_input('mode1_error_compare_fullFisher')
mode1_error_compare_heuristic = load_input('mode1_error_compare_heuristic')
mode1_MLE_compare_heuristic = load_input('mode1_MLE_compare_heuristic')
mode1_MLE_compare_numeric = load_input('mode1_MLE_compare_numeric')
mode1_MLE_compare_prior = load_input('mode1_MLE_compare_prior')
mode1_composite_Err_N_Black_Lower = load_input('mode1_composite_Err_N_Black_Lower')
mode1_composite_Err_N_Black_Upper = load_input('mode1_composite_Err_N_Black_Upper')
mode1_composite_Err_N_Red_Lower = load_input('mode1_composite_Err_N_Red_Lower')
mode1_composite_Err_N_Red_Upper = load_input('mode1_composite_Err_N_Red_Upper')
mode1_composite_mean_N_Black = load_input('mode1_composite_mean_N_Black')
mode1_composite_mean_N_Red = load_input('mode1_composite_mean_N_Red')
