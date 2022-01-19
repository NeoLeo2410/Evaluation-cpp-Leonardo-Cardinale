#include "Matrix.cpp"

unsigned nx_Q2 = 21; // Nombre de points

double dx_Q2 = 1.0/(nx_Q2-1); // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector_Q2(double& dx){
    std::vector<double> initial(nx_Q2,0.0);
    for (unsigned i = 0; i < nx_Q2; i++){
        double x = i * dx;
        double y;
        y = 0.5 + std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial[i] = y;
    }
    return initial;
}

std::vector<double> vec_Q2 = initial_vector_Q2(dx_Q2);
//unsigned n_Q2 = vec_Q2.size();
Matrix D_Q2(1,nx_Q2,1.0); // Matrice de conductivité thermique

Matrix K_Q2(Matrix D, double dx){
    unsigned n = D_Q2.getcolumns();
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

Matrix f_Q2(Matrix T){ 
    Matrix T1(T.transpose());
    Matrix m(K_Q2(D_Q2,dx_Q2) * T1);
    Matrix n(m.transpose());
    return n;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_explicite_Q2(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec_Q2};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f_Q2(prev) * step);
        solution.push_back((prev + next).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_Q2(std::vector<std::vector<double> > v){
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

/* int main(){
    double horiz = 0.5; // Horizon temporelle
    double dt = horiz/(nt-1); // Pas temporel
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution = euler_explicite(dt,horiz);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << " μs" << std::endl;
    exportsolution(solution);
    K(D,dx).print();
    return EXIT_SUCCESS;
} */
