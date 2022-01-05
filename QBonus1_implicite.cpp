#include "Matrix.cpp"
#include <fstream>

double dx = 0.01; // Pas spatial

// Construction de la condition initiale

std::vector<double> initial_vector(double& dx){
    std::vector<double> initial;
    initial.push_back(0.0);
    double x = 0;
    while (x + dx <= 1){
        x += dx;
        double y;
        y = 0.5 + 0.5 * std::sin(2 * M_PI * x) - 0.5 * std::cos(2 * M_PI * x);
        initial.push_back(y);
    }
    return initial;
}

std::vector<double> vec = initial_vector(dx);
unsigned n = vec.size();

// On travaille ici avec une matrice de conductivité aléatoire, de valeurs comprises entre 0.5 et 1.5 (en unités SI)

Matrix cond(){
    std::vector<double> dvec;
    double r;
    for (unsigned i = 0; i < n; i++){
        r = 0.5 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(1.5-0.5)));
        dvec.push_back(r);
    }
    Matrix d(dvec,0);
    return d;
}

Matrix D(cond()); // Matrice de conductivité thermique

Matrix K(Matrix D, double dx){
    unsigned n = D.getcolumns();
    Matrix cond(n, n, 0.0);
    for (unsigned i = 0; i < n; i++){
        if (i+1 < n){
            cond(i,i+1) = D(0,i+1)/(std::pow(dx,2));
            cond(i,i) = - (D(0,i) + D(0,i+1))/std::pow(dx,2);
        }
        else if (i-1 >= 0){
            cond(i,i-1) = D(0,i)/std::pow(dx,2);
        }
    }
    return cond;
}

// Fonction qui entre en jeu dans la méthode d'Euler. On travaille avec des vecteurs convertis en matrices lignes, d'où le passage par la transposée pour la phase de calcul, puis de nouveau pour retourner le résultat

Matrix f(Matrix T){ 
    Matrix T1(T.transpose());
    Matrix m((K(D,dx) * T1) * (-1.0));
    Matrix n(m.transpose());
    return n;
}

// Pivot de Gauss. Pas convaincu que ce soit la bonne approche pour l'instant pour des raisons de complexité

Matrix gauss(Matrix A, Matrix B){
    unsigned h = 0;
    unsigned k = 0;
    unsigned n = A.getrows();
    while ((h <  n) && (k < n)){
        unsigned i_max = h;
        for (unsigned l = h+1; l < n; l++){
            if (std::abs(A(l,k)) > A(i_max,k)){
                i_max = l;
            }
        }
        if (A(i_max,k) == 0){
            k++;
        }
        else{
            A.swap(h,i_max);
            for (unsigned i = h+1; i < n; i++){
                double f = A(i,k)/A(h,k);
                A(i,k) = 0;
                for (unsigned j = k+1; j < n; j++){
                    A(i,j) = A(i,j) - A(h,j) * f;
                    B(0,i) = B(0,i) - A(h,j) * f;
                }
            }
            h++;
            k++;
        }
    }
    Matrix X(1,n,0.0);
    X(0,n-1) = B(0,n-1)/A(n-1,n-1);
    unsigned m = n-2;
    while (m >= 0){
        double x = B(0,m);
        for (unsigned p = m+1; p < n; p++){
            x -= A(m,p) * X(0,p);
        }
        X(0,m) = x/A(m,m);
    }
    return X;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_implicite(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f(prev) * step);
        Matrix x0(prev + next);
        Matrix I(n,n,1.0);
        solution.push_back(gauss(K(D,dx) - I,prev*(-1.0)).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    std::ofstream myfile;
    myfile.open("QBonus1_implicite.txt");
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
    double dt = 0.01; // Pas temporel
    double horiz = 0.5; // Horizon temporelle
    std::vector<std::vector<double> > solution = euler_implicite(dt,horiz);
    exportsolution(solution);
    return EXIT_SUCCESS;
}