import csv
import pandas
from bokeh.plotting import figure
from bokeh.io import output_file, show
import glob
'''
Created by Sebastian Pajak
Solve set of equations
-2x1 + x2 = 0
x1 + -2x2 +x3  = 0
x2 +-2x3 +x4 =
1x3 - 2x4 = 1
'''

def define_matrix():
    matrix = [[-2, 1, 0, 0, 0],
              [1, -2, 1, 0, 0],
              [0, 1, -2, 1, 0],
              [0, 0, 1, -2, -1]]

    return matrix


def jacobi_method(matrix):
    i = 0
    x1_old = 0
    x2_old = 0
    x3_old = 0
    x4_old = 0
    x1_new = x1_old
    x2_new = x2_old
    x3_new = x3_old
    x4_new = x4_old
    res = 1
    residua_to_plot = {}

    while res > 0.0000001:
        residua_to_plot.update({i: res})
        x1_new = (matrix[0][4] - matrix[0][1]*x2_old - matrix[0][2]*x3_old - matrix[0][3]*x4_old) / matrix[0][0]
        x2_new = (matrix[1][4] - matrix[1][0]*x1_old - matrix[1][2]*x3_old - matrix[1][3]*x4_old) / matrix[1][1]
        x3_new = (matrix[2][4] - matrix[2][0]*x1_old - matrix[2][1]*x2_old - matrix[2][3]*x4_old) / matrix[2][2]
        x4_new = (matrix[3][4] - matrix[3][0]*x1_old - matrix[3][1]*x2_old - matrix[3][2]*x3_old) / matrix[3][3]
        res = abs(x1_new + x2_new + x3_new + x4_new - x1_old - x2_old - x3_old - x4_old)

        x1_old = x1_new
        x2_old = x2_new
        x3_old = x3_new
        x4_old = x4_new

        i += 1

    with open('data.csv', 'a+', newline='') as csvfile:
        field_names = ['Method', 'Iterations']
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writeheader()
        writer.writerow({'Method': 'Jacobi', 'Iterations': str(i)})
    with open('jacobi_residua.csv','w', newline='') as csvfile:
        field_names = ['Iteration', 'Residuum']
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writeheader()
        for k,v in residua_to_plot.items():
            writer.writerow({'Iteration': k, 'Residuum':v})

    return(x1_new, x2_new, x3_new, x4_new)

def gauss_siedl(matrix):
    i = 0
    x1_old = 0
    x2_old = 0
    x3_old = 0
    x4_old = 0
    x1_new = x1_old
    x2_new = x2_old
    x3_new = x3_old
    x4_new = x4_old
    res = 1
    residua_to_plot = {}

    while res > 0.0000001:
        residua_to_plot.update({i: res})
        x1_new = (matrix[0][4] - matrix[0][1]*x2_old - matrix[0][2]*x3_old - matrix[0][3]*x4_old) / matrix[0][0]
        x2_new = (matrix[1][4] - matrix[1][0]*x1_new - matrix[1][2]*x3_old - matrix[1][3]*x4_old) / matrix[1][1]
        x3_new = (matrix[2][4] - matrix[2][0]*x1_new - matrix[2][1]*x2_new - matrix[2][3]*x4_old) / matrix[2][2]
        x4_new = (matrix[3][4] - matrix[3][0]*x1_new - matrix[3][1]*x2_new - matrix[3][2]*x3_new) / matrix[3][3]
        res = abs(x1_new + x2_new + x3_new + x4_new - x1_old - x2_old - x3_old - x4_old)

        x1_old = x1_new
        x2_old = x2_new
        x3_old = x3_new
        x4_old = x4_new

        i += 1
    with open('data.csv', 'a+', newline='') as csvfile:
        field_names = ['Method', 'Iterations']
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writerow({'Method': 'Gauss-Siedl', 'Iterations': str(i)})
    with open('gauss_siedl_residua.csv','w', newline='') as csvfile:
        field_names = ['Iteration', 'Residuum']
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writeheader()
        for k,v in residua_to_plot.items():
            writer.writerow({'Iteration': k, 'Residuum':v})


    return(x1_new, x2_new, x3_new, x4_new)

def sor(matrix):
    i = 0
    x1_old = 0
    x2_old = 0
    x3_old = 0
    x4_old = 0
    x1_new = x1_old
    x2_new = x2_old
    x3_new = x3_old
    x4_new = x4_old
    res = 1
    omega = 1.2
    residua_to_plot = {}

    while res > 0.0000001:
        residua_to_plot.update({i: res})
        x1_new = omega*(matrix[0][4] - matrix[0][1]*x2_old - matrix[0][2]*x3_old - matrix[0][3]*x4_old) / matrix[0][0] + (1-omega)*x1_old
        x2_new = omega*(matrix[1][4] - matrix[1][0]*x1_new - matrix[1][2]*x3_old - matrix[1][3]*x4_old) / matrix[1][1] + (1-omega)*x2_old
        x3_new = omega*(matrix[2][4] - matrix[2][0]*x1_new - matrix[2][1]*x2_new - matrix[2][3]*x4_old) / matrix[2][2] + (1-omega)*x3_old
        x4_new = omega*(matrix[3][4] - matrix[3][0]*x1_new - matrix[3][1]*x2_new - matrix[3][2]*x3_new) / matrix[3][3] + (1-omega)*x4_old
        res = abs(x1_new + x2_new + x3_new + x4_new - x1_old - x2_old - x3_old - x4_old)

        x1_old = x1_new
        x2_old = x2_new
        x3_old = x3_new
        x4_old = x4_new

        i += 1

    with open('data.csv', 'a+', newline='') as csvfile:
        field_names = ['Method', 'Iterations']
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writerow({'Method': 'SOR 1.2', 'Iterations': str(i)})
    with open('sor_residua.csv','w', newline='') as csvfile:
        field_names = ['Iteration', 'Residuum']
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writeheader()
        for k,v in residua_to_plot.items():
            writer.writerow({'Iteration': k, 'Residuum':v})

    return(x1_new, x2_new, x3_new, x4_new)

def plotting_data():
    #add to list all files with extensions .csv
    csv_data = []
    for file in glob.glob("*.csv"):
        csv_data.append(file)

    #import data from files
    df_res_jacobi = pandas.read_csv(csv_data[2])
    df_res_gauss_siedl = pandas.read_csv(csv_data[1])
    df_res_sor = pandas.read_csv(csv_data[3])
    df_compare_methods = pandas.read_csv(csv_data[0])

    #data for plotting
    x1 = df_res_jacobi["Iteration"]
    y1 = df_res_jacobi["Residuum"]

    x2 = df_res_gauss_siedl["Iteration"]
    y2 = df_res_gauss_siedl["Residuum"]

    x3 = df_res_sor["Iteration"]
    y3 = df_res_sor["Residuum"]

    x4 = df_compare_methods["Method"]
    y4 = df_compare_methods["Iterations"]

    # output to static HTML file
    output_file("iteration_methods.html", title="line plot example")
    p4 = figure(x_range=x4, plot_height=350, title="Iterations per method")
    p1 = figure(width=500, height=300, title="Jacobi")
    p2 = figure(width=500, height=300, title="Gauss-Siedl")
    p3 = figure(width=500, height=300, title="SOR 1.2")
    p1.line(x1, y1)
    p1.xaxis.axis_label = "Iterations"
    p1.yaxis.axis_label = "Residuum"
    p2.line(x2, y2, color="red")
    p2.xaxis.axis_label = "Iterations"
    p2.yaxis.axis_label = "Residuum"
    p3.line(x3, y3, color="green")
    p3.xaxis.axis_label = "Iterations"
    p3.yaxis.axis_label = "Residuum"
    p4.vbar(x=x4, top=y4, width=0.05, color="Black")
    p4.xgrid.grid_line_color = None
    show(p1)
    show(p2)
    show(p3)
    show(p4)


if __name__ == "__main__":
    matrix = define_matrix()
    jacobi_method(matrix)
    gauss_siedl(matrix)
    sor(matrix)
    plotting_data()