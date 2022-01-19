#include "QBonus1_implicite.cpp"
#include "Eigen/Sparse"

unsigned nx_QB2e = 21; // Nombre de points

double dx_QB2e = 1.0/(nx_QB2e-1); // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector_QB2e(double& dx){
    std::vector<double> initial(nx_QB2e,0.0);
    for (unsigned i = 0; i < nx_QB2e; i++){
        double x = i * dx;
        double y;
        y = 0.5 + std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial[i] = y;
    }
    return initial;
}

std::vector<double> vec_QB2e = initial_vector_QB2e(dx_QB2e);
Matrix D_QB2e(1,nx_QB2e,1.0); // Matrice de conductivité thermique

Eigen::SparseMatrix<double> K_QB2e(Matrix D, double dx){
    unsigned n = D.getcolumns();
    Eigen::SparseMatrix<double> cond(n, n);
    for (unsigned i = 0; i < n; i++){
        if (i+1 < n){
            cond.coeffRef(i,i+1) = D(0,i+1)/(std::pow(dx,2));
            cond.coeffRef(i,i) = - (D(0,i) + D(0,i+1))/std::pow(dx,2);
        }
        if (i != 0){
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

Eigen::VectorXd f_QB2e(Eigen::VectorXd T){
    return K_QB2e(D_QB2e,dx_QB2e) * T;
}

// Conversion Matrix -> Eigen::MatrixXd

Eigen::MatrixXd matrixtoeigen_QB2e(Matrix M){
    unsigned n = M.getrows();
    unsigned p = M.getcolumns();
    Eigen::MatrixXd m(n,p);
    for (unsigned i = 0; i < n; i++){
        for (unsigned j = 0; j < p; j++){
            m(i,j) = M(i,j);
        }
    } 
    return m;
}

// Conversion std::vector -> Eigen::VectorXd

Eigen::VectorXd vectortoeigen_QB2e(std::vector<double> vec){
    unsigned n = vec.size();
    Eigen::VectorXd V(n);
    for (unsigned i = 0; i < n; i++){
        V(i) = vec[i];
    }
    return V;
}

// Conversion Eigen::VectorXd -> std::vector

std::vector<double> eigentovector_QB2e(Eigen::VectorXd v){
    unsigned n = v.rows();
    std::vector<double> vec(n);
    for (unsigned i = 0; i < n; i++){
        vec[i] = v(i);
    }
    return vec;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_explicite_QB2e(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec_QB2e};
    while (dates[dates.size() - 1] + step < T){
        Eigen::VectorXd prev(nx_QB2e);
        prev = vectortoeigen_QB2e(solution[solution.size() - 1]);
        Eigen::VectorXd next(nx_QB2e);
        next = f_QB2e(prev) * step;
        solution.push_back(eigentovector_QB2e(prev + next));
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_QB2e(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    unsigned m = v[0].size();
    std::ofstream myfile;
    myfile.open("QBonus2_explicite.txt");
    for (unsigned i = 0; i < n - 1; i++){
        for (unsigned j = 0; j < m - 1; j++){
            myfile << v[i][j] << " ";
        }
        myfile << v[i][m-1] << "\n";
    }
    for (unsigned j = 0; j < m - 1; j++){
        myfile << v[n-1][j] << " ";
    }
    myfile << v[n-1][m-1] << std::endl;
    myfile.close();
}

/* int main(){
    double horiz = 0.5; // Horizon temporelle
    double dt = horiz/(nt-1); // Pas temporel
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution = euler_explicite(dt,horiz);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << " μs" << std::endl;
    exportsolution(solution);
    std::cout << Eigen::MatrixXd(K(D,dx)) << std::endl;
    return EXIT_SUCCESS;
} */
