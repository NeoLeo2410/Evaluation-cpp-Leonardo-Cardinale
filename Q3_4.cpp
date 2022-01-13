#include "Matrix.cpp"
#include <fstream>

unsigned nx = 21; // Nombre de points
unsigned nt = 1001; // Nombre de dates

double dx = 1.0/nx; // Pas spatial

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

Matrix f(Matrix T){ 
    Matrix T1(T.transpose());
    Matrix m(K(D,dx) * T1);
    Matrix n(m.transpose());
    return n;
}

// On utilise la norme (au carré) dans l'algorithme du gradient conjugué qui suit

double norm_squared(Matrix& A){
    double result = 0.0;
    unsigned n = A.getrows();
    unsigned i;
    unsigned j;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            result += std::pow(std::abs(A(i,j)),2);
        }
    }
    return result;
}

// Pour résoudre des systèmes linéaires associés à une matrice symétrique définie positive

Matrix gradient_conjugue(Matrix A, Matrix b, Matrix x0, double eps){
    Matrix X(x0.transpose());
    Matrix r(b.transpose() - (A * X));
    Matrix p(r);
    unsigned k = 0;
    unsigned n = X.getrows();
    double rr = (r.transpose() * r)(0,0);
    while ((norm_squared(r) > std::pow(eps,2)) && (k < n)){
        Matrix Ap(A * p);
        double denom = (p.transpose() * Ap)(0,0);
        double alpha = rr/denom;
        X = X + p * alpha;
        r = r - Ap * alpha;
        double rrtemp = (r.transpose() * r)(0,0);
        double beta = rrtemp/rr;
        rr = rrtemp;
        p = r + p * beta;
        k++;
    }
    return X.transpose();
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_implicite(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f(prev) * step);
        Matrix x0(prev + next);
        Matrix I(n,n,0.0); // construction sur-le-champ d'une matrice identité
        for (unsigned i = 0; i < n; i++){
            I(i,i) = 1.0;
        }
        solution.push_back(gradient_conjugue(I - (K(D,dx)*step),prev,x0,dx).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    std::ofstream myfile;
    myfile.open("Q3_4.txt");
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
    double dt = horiz/nt; // Pas temporel
    std::vector<std::vector<double> > solution = euler_implicite(dt,horiz);
    exportsolution(solution);
    return EXIT_SUCCESS;
}