#include "QBonus2_explicite.cpp"
#include "Eigen/SparseLU"

unsigned nx_QB2i = 21; // Nombre de points

double dx_QB2i = 1.0/(nx_QB2i-1); // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector_QB2i(double& dx){
    std::vector<double> initial(nx_QB2i,0.0);
    for (unsigned i = 0; i < nx_QB2i; i++){
        double x = i * dx;
        double y;
        y = 0.5 + std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial[i] = y;
    }
    return initial;
}

std::vector<double> vec_QB2i = initial_vector_QB2i(dx_QB2i);
Matrix D_QB2i(1,nx_QB2i,1.0); // Matrice de conductivité thermique

Eigen::SparseMatrix<double> K_QB2i(Matrix D, double dx){
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

Eigen::VectorXd f_QB2i(Eigen::VectorXd T){
    return K_QB2i(D_QB2i,dx_QB2i) * T;
}

// Conversion Matrix -> Eigen::MatrixXd

Eigen::MatrixXd matrixtoeigen_QB2i(Matrix M){
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

Eigen::VectorXd vectortoeigen_QB2i(std::vector<double> vec){
    unsigned n = vec.size();
    Eigen::VectorXd V(n);
    for (unsigned i = 0; i < n; i++){
        V(i) = vec[i];
    }
    return V;
}

// Conversion Eigen::VectorXd -> std::vector

std::vector<double> eigentovector_QB2i(Eigen::VectorXd v){
    unsigned n = v.rows();
    std::vector<double> vec(n);
    for (unsigned i = 0; i < n; i++){
        vec[i] = v(i);
    }
    return vec;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_implicite_QB2i(const double& step, const double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec_QB2i};
    while (dates[dates.size() - 1] + step < T){
        Eigen::VectorXd prev(nx_QB2i);
        prev = vectortoeigen_QB2i(solution[solution.size() - 1]);
        Eigen::SparseMatrix<double> I(nx_QB2i,nx_QB2i);
        for (unsigned i = 0; i < nx_QB2i; i++){
            I.coeffRef(i,i) = 1.0;
        }
        Eigen::SparseMatrix<double> A;
        A = I - (K_QB2i(D_QB2i,dx_QB2i)*step);
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        solution.push_back(eigentovector_QB2i(solver.solve(prev)));
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_QB2i(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    unsigned m = v[0].size();
    std::ofstream myfile;
    myfile.open("QBonus2_implicite.txt");
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