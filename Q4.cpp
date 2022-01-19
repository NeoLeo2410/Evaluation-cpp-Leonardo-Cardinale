#include "Q2.cpp"

// On utilise la norme (au carré) dans l'algorithme du gradient conjugué qui suit.

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

std::vector<std::vector<double> > euler_implicite_Q4(double& step, double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f(prev) * step);
        Matrix x0(prev + next);
        Matrix I(nx,nx,0.0); // construction sur-le-champ d'une matrice identité
        for (unsigned i = 0; i < nx; i++){
            I(i,i) = 1.0;
        }
        solution.push_back(gradient_conjugue(I - (K(D,dx)*step),prev,x0,dx).tovector());
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_Q4(std::vector<std::vector<double> > v){
    unsigned n = v.size();
    unsigned m = v[0].size();
    std::ofstream myfile;
    myfile.open("Q4.txt");
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