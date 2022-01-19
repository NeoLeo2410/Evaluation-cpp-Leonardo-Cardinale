#include "QBonus1_explicite.cpp"
#include "Eigen/Dense"

// Conversion Matrix -> Eigen::MatrixXd

Eigen::MatrixXd matrixtoeigen(Matrix M){
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

// Conversion matrice ligne -> Eigen::VectorXd

Eigen::VectorXd vectortoeigen(Matrix V){
    unsigned n = V.getrows();
    unsigned p = V.getcolumns();
    if (n != 1){
        throw std::runtime_error("Argument must be row matrix");
    }
    else{
        Eigen::VectorXd v(p);
        for (unsigned j = 0; j < p; j++){
            v(j) = V(0,j);
        }
    return v;
    }
}

// Conversion Eigen::VectorXd -> std::vector

std::vector<double> eigentovector(Eigen::VectorXd v){
    unsigned n = v.rows();
    std::vector<double> vec(n);
    for (unsigned i = 0; i < n; i++){
        vec[i] = v(i);
    }
    return vec;
}

// Avec la méthode d'Euler, on obtient un vecteur de vecteurs où l'élément (i,j) représente T_j(i*dt)

std::vector<std::vector<double> > euler_implicite_QB1i(const double& step, const double& T){
    std::vector<double> dates {0.0};
    std::vector<std::vector<double> > solution {vec};
    while (dates[dates.size() - 1] + step < T){
        Matrix prev(solution[solution.size() - 1],0);
        Matrix next(f_rand(prev) * step);
        Matrix x0(prev + next);
        Matrix I(nx,nx,0.0);
        for (unsigned i = 0; i < nx; i++){
            I(i,i) = 1.0;
        }
        Eigen::MatrixXd A(nx,nx);
        A = matrixtoeigen(I - (K(D_rand,dx)*step));
        Eigen::VectorXd B(nx);
        B = vectortoeigen(prev);
        solution.push_back(eigentovector(A.lu().solve(B)));
        dates.push_back(dates[dates.size() - 1] + step);
    }
    return solution;
}

// Pour exporter au format .txt une liste de listes pouvant être passée en argument à numpy.array() en Python

void exportsolution_QB1i(std::vector<std::vector<double> >& v){
    unsigned n = v.size();
    unsigned m = v[0].size();
    std::ofstream myfile;
    myfile.open("QBonus1_implicite.txt");
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