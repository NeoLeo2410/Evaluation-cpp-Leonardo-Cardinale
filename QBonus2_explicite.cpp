#include "QBonus1_implicite.cpp"
#include "Eigen/Sparse"

// Matrice de diffusion

Eigen::SparseMatrix<double> K_sparse(Matrix D, double dx){
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

Eigen::VectorXd f_sparse(Eigen::VectorXd T){
    return K_sparse(D,dx) * T;
}

// Conversion std::vector -> Eigen::VectorXd

Eigen::VectorXd vectortoeigen1(std::vector<double> vec){
    unsigned n = vec.size();
    Eigen::VectorXd V(n);
    for (unsigned i = 0; i < n; i++){
        V(i) = vec[i];
    }
    return V;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_explicite_QB2e(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Eigen::VectorXd prev(nx);
        prev = vectortoeigen1(solution[solution.size() - 1]);
        Eigen::VectorXd next(nx);
        next = f_sparse(prev) * step;
        solution.push_back(eigentovector(prev + next));
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