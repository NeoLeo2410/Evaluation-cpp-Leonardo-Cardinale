#include "Matrix.cpp"

const unsigned nx = 21; // Nombre de points

const double dx = 1.0/(nx-1); // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector(const double& dx){
    std::vector<double> initial(nx,0.0);
    for (unsigned i = 0; i < nx; i++){
        double x = i * dx;
        double y;
        y = 0.5 + std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial[i] = y;
    }
    return initial;
}

const std::vector<double> vec = initial_vector(dx);
Matrix D(1,nx,1.0); // Matrice de conductivité thermique

Matrix K(Matrix& D, const double& dx){
    unsigned n = D.getcolumns();
    Matrix cond(n, n, 0.0);
    for (unsigned i = 0; i < n; i++){
        if (i+1 < n){
            cond(i,i+1) = D(0,i+1)/(std::pow(dx,2));
            cond(i,i) = - (D(0,i) + D(0,i+1))/std::pow(dx,2);
        }
        if (i != 0){
            cond(i,i-1) = D(0,i)/std::pow(dx,2);
        }
    }
    for (unsigned i = 1; i < n; i++){
        cond(0,i) = 0;
    }
    for (unsigned i = 0; i < n - 1; i++){
        cond(n-1,i) = 0;
    }
    cond(0,0) = - (D(0,0) + D(0,1))/std::pow(dx,2);
    cond(n-1,n-1) = - (D(0,n-1) + D(0,n-1))/std::pow(dx,2);
    return cond;
}

// Fonction qui entre en jeu dans la méthode d'Euler. On travaille avec des vecteurs convertis en matrices lignes, d'où le passage par la transposée pour la phase de calcul, puis de nouveau pour retourner le résultat

Matrix f(Matrix& T){ 
    Matrix T1(T.transpose());
    Matrix m(K(D,dx) * T1);
    Matrix n(m.transpose());
    return n;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_explicite_Q2(const double& step, const double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f(prev) * step);
        solution.push_back((prev + next).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_Q2(std::vector<std::vector<double> >& v){
    unsigned n = v.size();
    unsigned m = v[0].size();
    std::ofstream myfile;
    myfile.open("Q2.txt");
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