#include "Matrix.cpp"

int main(){
    Matrix m(1,2,0.0);
    for (unsigned i = 0; i < 1; i++){
        for (unsigned j = 0; j < 2; j++){
            if (i == 0){
                m(i,j) = j + 1;
            }
            else{
                m(i,j) = 3 + j;
            }
        }
    }
    Matrix n(1,2,5.0);
    Matrix l(m - n);
    m.print();
    n.print();
    l.print();
    std::vector<double> vec {0.5, 1.0, 4.0};
    vec.push_back(0.0);
    vec.push_back(42.0);
    vec.push_back(69.0);
    Matrix M(vec,0);
    M.print();
    Matrix N(M.transpose());
    N.print();
    Matrix P(1,10,0.231);
    std::vector<double> vect = P.tovector();
    std::cout << vect[0] << std::endl;
    return EXIT_SUCCESS;
}