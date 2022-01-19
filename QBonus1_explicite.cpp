#include "Q3_4.cpp"
#include <random>

unsigned nx_QB1e = 21; // Nombre de points

double dx_QB1e = 1.0/(nx_QB1e-1); // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector_QB1e(double& dx){
    std::vector<double> initial(nx_QB1e,0.0);
    for (unsigned i = 0; i < nx_QB1e; i++){
        double x = i * dx;
        double y;
        y = 0.5 + std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial[i] = y;
    }
    return initial;
}

std::vector<double> vec_QB1e = initial_vector_QB1e(dx_QB1e);

// On travaille ici avec une matrice de conductivité aléatoire, de valeurs comprises entre 0.5 et 1.5 (en unités SI)

Matrix cond_QB1e(){
    std::vector<double> dvec;
    double r;
    for (unsigned i = 0; i < nx_QB1e; i++){
        r = 0.5 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(1.5-0.5)));
        dvec.push_back(r);
    }
    Matrix d(dvec,0);
    return d;
}

Matrix D_QB1e(cond_QB1e()); // Matrice de conductivité thermique

Matrix K_QB1e(Matrix D, double dx){
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

Matrix f_QB1e(Matrix T){ 
    Matrix T1(T.transpose());
    Matrix m(K_QB1e(D_QB1e,dx_QB1e) * T1);
    Matrix n(m.transpose());
    return n;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_explicite_QB1e(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec_QB1e};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f_QB1e(prev) * step);
        solution.push_back((prev + next).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_QB1e(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    unsigned m = v[0].size();
    std::ofstream myfile;
    myfile.open("QBonus1_explicite.txt");
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
    std::vector<std::vector<double> > solution = euler_explicite(dt,horiz);
    exportsolution(solution);
    K(D,dx).print();
    return EXIT_SUCCESS;
} */
