#include "Matrix.cpp"
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

unsigned nx = 21; // Nombre de points
unsigned nt = 1001; // Nombre de dates

double dx = 1.0/(nx-1); // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector(double& dx){
    std::vector<double> initial(nx,0.0);
    for (unsigned i = 0; i < nx; i++){
        double x = i * dx;
        double y;
        y = 0.5 + 0.5 * std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial[i] = y;
    }
    return initial;
}

std::vector<double> vec = initial_vector(dx);
unsigned n = vec.size();
Matrix D(1,n,1.0); // Matrice de conductivité thermique

Eigen::SparseMatrix<double> K(Matrix D, double dx){
    unsigned n = D.getcolumns();
    Eigen::SparseMatrix<double> cond(n, n);
    for (unsigned i = 0; i < n; i++){
        if (i+1 < n){
            cond.coeffRef(i,i+1) = D(0,i+1)/(std::pow(dx,2));
            cond.coeffRef(i,i) = - (D(0,i) + D(0,i+1))/std::pow(dx,2);
        }
        else if (i-1 >= 0){
            cond.coeffRef(i,i-1) = D(0,i)/std::pow(dx,2);
        }
    }
    for (unsigned i = 1; i < n; i++){
        cond.coeffRef(0,i) = 0;
    }
    for (unsigned i = 0; i < n - 1; i++){
        cond.coeffRef(n-1,i) = 0;
    }
    cond.coeffRef(0,0) = - (D(0,0) + D(0,1))/std::pow(dx,2);
    cond.coeffRef(n-1,n-1) = - (D(0,n-1) + D(0,n-1))/std::pow(dx,2);
    return cond;
}

// Fonction qui entre en jeu dans la méthode d'Euler. On travaille avec des vecteurs convertis en matrices lignes, d'où le passage par la transposée pour la phase de calcul, puis de nouveau pour retourner le résultat

Eigen::VectorXd f(Eigen::VectorXd T){
    return K(D,dx) * T;
}

// Conversion Eigen::VectorXd -> std::vector

std::vector<double> eigentovector(Eigen::VectorXd v){
    unsigned n = v.rows();
    std::vector<double> vec(n);
    for (unsigned i = 0; i < n; i++){
        vec[i] = v(i);
    }
    return vec;
}

// Conversion std::vector -> Eigen::VectorXd

Eigen::VectorXd vectortoeigen(std::vector<double> vec){
    unsigned n = vec.size();
    Eigen::VectorXd V(n);
    for (unsigned i = 0; i < n; i++){
        V(i) = vec[i];
    }
    return V;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_implicite(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Eigen::VectorXd prev(n);
        prev = vectortoeigen(solution[solution.size() - 1]);
        Eigen::SparseMatrix<double> I(n,n);
        for (unsigned i = 0; i < n; i++){
            I.coeffRef(i,i) = 1.0;
        }
        Eigen::SparseMatrix<double> A;
        A = I - (K(D,dx)*step);
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
        solution.push_back(eigentovector(solver.compute(A).solve(prev)));
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    std::ofstream myfile;
    myfile.open("QBonus2_implicite.txt");
    myfile << "[";
    for (unsigned i = 0; i < n - 1; i++){
        myfile << "[";
        for (unsigned j = 0; j < v[i].size() - 1; j++){
            myfile << v[i][j] << ",";
        }
        myfile << v[i][v[i].size() - 1] << "],";
    }
    myfile << "[";
    for (unsigned j = 0; j < v[n-1].size() - 1; j++){
        myfile << v[n-1][j] << ",";
    }
    myfile << v[n-1][v[n-1].size() - 1] << "]]" << std::endl;
    myfile.close();
}

int main(){
    double horiz = 0.5; // Horizon temporelle
    double dt = horiz/(nt-1); // Pas temporel
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution = euler_implicite(dt,horiz);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << " μs" << std::endl;
    exportsolution(solution);
    return EXIT_SUCCESS;
}