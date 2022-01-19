#include "Q4.cpp"
#include <random>

// On travaille ici avec une matrice de conductivité aléatoire, de valeurs comprises entre 0.5 et 1.5 (en unités SI)

Matrix cond(){
    std::vector<double> dvec;
    double r;
    for (unsigned i = 0; i < nx; i++){
        r = 0.5 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(1.5-0.5)));
        dvec.push_back(r);
    }
    Matrix d(dvec,0);
    return d;
}

Matrix D_rand(cond()); // Matrice de conductivité thermique

// Fonction qui entre en jeu dans la méthode d'Euler. On travaille avec des vecteurs convertis en matrices lignes, d'où le passage par la transposée pour la phase de calcul, puis de nouveau pour retourner le résultat

Matrix f_rand(Matrix& T){ 
    Matrix T1(T.transpose());
    Matrix m(K(D_rand,dx) * T1);
    Matrix n(m.transpose());
    return n;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_explicite_QB1e(const double& step, const double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f_rand(prev) * step);
        solution.push_back((prev + next).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_QB1e(std::vector<std::vector<double> >& v){
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