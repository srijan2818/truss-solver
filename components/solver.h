#pragma once
#include <iostream>
#include <vector>  // for creating matrices nd vectors
#include <cmath>
#include <string>
#include <iomanip> // for setw formatting in printing
#include <algorithm> // for swapping rows only
#if defined(_MSC_VER)
    #define _RESTRICT __restrict
#else
    #define _RESTRICT __restrict__
#endif

class Matrix{

private:
    int rows,cols;
    std::vector<double> data;

public:
    Matrix(int r,int c): rows(r),cols(c),data(r*c,0.0) {}
    
    inline double& operator()(int i,int j){
        return data[i*cols+j];
    };
    inline const double& operator()(int i,int j) const{
        return data[i*cols+j];
    };

    int nrows() const{
        return rows;
    };
    int ncols() const{
        return cols;
    };
    double* raw(){
        return data.data();
    };
    const double* raw() const {
        return data.data();
    };

    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
    for (int i = 0; i < mat.nrows(); ++i) {
        os << "[ ";
        for (int j = 0; j < mat.ncols(); ++j) {
            os << std::setw(8) << mat(i, j) << " ";
        }
        os << "]\n";
    }
    return os;
}
};

inline double dot(double* _RESTRICT x,double* _RESTRICT y,int n){
    double s=0;
    for(int i=0;i<n;i++){
        s += x[i]*y[i];
    }
    return s;
}
inline void scalar_mult_add(double alpha,double* _RESTRICT x,double* _RESTRICT y,int n){
    for(int i=0;i<n;i++){
        y[i] += alpha*x[i];
    }
}
inline void matvec(const Matrix& A,double* _RESTRICT x,double* _RESTRICT y){
    const double* ptr = A.raw();
    const int rows = A.nrows();
    const int cols = A.ncols();
    for(int i=0;i<rows;i++){
        double sum= 0.0;
        const double* row_iset = ptr + i*cols;
        for(int j=0;j<cols;j++){
            sum += row_iset[j] * x[j];
        }
        y[i] = sum;
    } 
}

inline void scalar_add(int n,double alpha,double* _RESTRICT x){
    for(int i=0;i<n;i++){
        x[i] += alpha;
    }
}

inline void forward_sub(
    const Matrix& L,
    const double* _RESTRICT b,
    double* _RESTRICT y
){
    const int n = L.nrows();
    const int cols = L.ncols();
    const double* ptr = L.raw();
    for(int i=0;i<n;i++){
        const double* row_i = ptr + i*cols;
        double sum = 0.0;
        if(std::abs(row_i[i])<1e-12) throw std::runtime_error("Zero on diagonal in forward sub"); 
        for(int j=0;j<i;j++){
            sum += row_i[j]*y[j];
        }
        y[i] = (b[i] - sum)/row_i[i];
    }
}


inline void back_sub(
    const Matrix&   U,
    const double* _RESTRICT b,
    double* _RESTRICT y
){
    const int n = U.nrows();
    const int cols = U.ncols();
    const double* ptr = U.raw();

    for(int i=n-1;i>=0;i--){
        double sum=0.0;
        const double* row_i = ptr + i*cols;
        if(std::abs(row_i[i])<1e-12) throw std::runtime_error("Zero on diagonal in backward sub.");
        for(int j=i+1;j<n;j++){
            sum += row_i[j] * y[j];
        }
        y[i] = (b[i] - sum)/row_i[i];
    }
}
inline void transpose(
    Matrix& A,
    Matrix& B
){
    const int n = A.nrows();
    const int cols = A.ncols();
    double* A_ptr = A.raw();
    double* B_ptr = B.raw();
    for(int i=0;i<n;i++){
        double* A_row_i = A_ptr + i*cols;
        for(int j=0;j<cols;j++){
            B_ptr[j*n+i] = A_row_i[j];
        }  
    }
}
// ~ O(2/3 n^3)
inline void gauss_eliimination_inplace(
    Matrix& A,
    double* b,
    double* x
){
    const int n = A.nrows();
    const int cols = A.ncols();
    double* A_ptr = A.raw();

    for(int i=0;i<n;i++){
        int pivot_row = i;
        double* A_row_i = A_ptr + i*cols;
        double max_val = std::abs(A_row_i[i]);
        for(int j=i+1;j<n;j++){
            if(std::abs(A_ptr[j*cols+i]) > max_val){
            max_val = std::abs(A_ptr[j*cols+i]);
            pivot_row=j;
            }
        }
        if(pivot_row != i){
            double* A_row_pivot = A_ptr + pivot_row*cols; 
            for(int k=0;k<cols;k++){
                std::swap(A_row_i[k],A_row_pivot[k]);
            }
            std::swap(b[i],b[pivot_row]);
        }
        const double pivot = A_row_i[i];
        for(int r=i+1;r<n;r++){
            double* A_row_r = A_ptr  + r*cols;
            const double factor = A_row_r[i]/pivot;
            for(int k=i;k<n;k++){
                A_row_r[k] -= factor * A_row_i[k];
            }
            b[r] -= factor*b[i];
        }
    }
    back_sub(A,b,x);
}
inline void gauss_elimination(
    const Matrix& A,
    const double* b,
    double* x
){
    Matrix A_work = A;
    const int n = A.nrows();
    std::vector<double> b_work(b, b + n);
    gauss_eliimination_inplace(A_work, b_work.data(), x);
}



// O(2/3 n^3)
inline void lu_doolittle(
    Matrix& A,
    Matrix& L,
    Matrix& U,
    std::vector<int>& P
){
    const int n = A.nrows();
    const int cols = A.ncols();
    P.resize(n);
    for(int i=0;i<n;i++) P[i] = i;
    double* A_ptr = A.raw();
    double* L_ptr = L.raw();
    double* U_ptr = U.raw();

    for(int i=0;i<n;i++){
        int pivot = i;
        const int row = i*cols;
        double* A_row_i = A_ptr + row;
        double* L_row_i = L_ptr + row;
        double* U_row_i = U_ptr + row;

        double max_val = std::abs(A_row_i[i]);
        for(int r=i+1;r<n;r++){
            if(std::abs(A_ptr[r*cols+i]) > max_val){
            max_val = std::abs(A_ptr[r*cols+i]);
            pivot =r;
            }
        }
        if(pivot != i){
            double* A_row_pivot = A_ptr + pivot*cols; 
            for(int k=0;k<cols;k++){
                std::swap(A_row_i[k],A_row_pivot[k]);
            }
            std::swap(P[i],P[pivot]);
            for(int k=0;k<i;k++){
            double* L_row_pivot = L_ptr + pivot*cols; 
                std::swap(L_row_i[k],L_row_pivot[k]);
            }

        }
        L_row_i[i] = 1.0;

        for(int j=i;j<n;j++){
            double sum = 0.0;
            for(int k=0;k<i;k++){
                sum += L_row_i[k] * U_ptr[k*cols+j];
            }
            U_row_i[j] = (A_row_i[j] - sum);
        }

        if(std::abs(U_row_i[i])<1e-12) {
            std::cout<<"Near 0 value at step / value "<<i<<U_row_i[i]<<std::endl;
            throw std::runtime_error("Singular / Near singular value.");}


        for(int j=i+1;j<n;j++){
            double sum=0.0;
            double* L_row_j = L_ptr + j*cols; 
            for(int k=0;k<i;k++){
                sum += L_row_j[k] * U_ptr[k*cols+i];
            }
            L_row_j[i] = (A_ptr[j*cols+i] - sum)/U_row_i[i];
        }

    }
}
// O(2/3 n^3)
inline void lu_crout(
    Matrix& A,
    Matrix& L,
    Matrix& U,
    std::vector<int>& P
){
    const int n = A.nrows();
    const int cols = A.ncols();
    P.resize(n);
    for(int i=0;i<n;i++) P[i] = i;
    double* A_ptr = A.raw();
    double* L_ptr = L.raw();
    double* U_ptr = U.raw();

    for(int i=0;i<n;i++){
        int pivot = i;
        const int row = i*cols;
        double* row_i = A_ptr + row;
        double* L_row_i = L_ptr + row;
        double* U_row_i = U_ptr + row;

        double max_val = std::abs(row_i[i]);
        for(int r=i+1;r<n;r++){
            if(std::abs(A_ptr[r*cols+i]) > max_val){
            max_val = std::abs(A_ptr[r*cols+i]);
            pivot =r;
            }
        }
        if(pivot != i){
            double* A_row_pivot = A_ptr + pivot*cols; 
            for(int k=0;k<cols;k++){
                std::swap(row_i[k],A_row_pivot[k]);
            }
            std::swap(P[i],P[pivot]);
            for(int k=0;k<i;k++){
            double* L_row_pivot = L_ptr + pivot*cols; 
                std::swap(L_row_i[k],L_row_pivot[k]);
            }

        }
        U_row_i[i] = 1.0;

        for(int j=i;j<n;j++){
            double sum = 0.0;
            double* L_row_j = L_ptr + j*cols;
            for(int k=0;k<i;k++){
                sum += L_row_j[k] * U_ptr[k*cols+i];
            }
            L_row_j[i] = (A_ptr[j*cols+i] - sum);
        }

        if(std::abs(L_row_i[i])<1e-12) {
            std::cout<<"Near 0 value at step / value "<<i<<L_row_i[i]<<std::endl;
            throw std::runtime_error("Singular / Near singular value.");}

        for(int j=i+1;j<n;j++){
            double sum=0.0;
            for(int k=0;k<i;k++){
                sum += L_row_i[k] * U_ptr[k*cols+j];
            }
            U_row_i[j] = (row_i[j] - sum)/L_row_i[i];
        }

    }
}
inline void lu_cholesky(
    Matrix& A,
    Matrix& L
){
    const int n = A.nrows();
    const int cols = A.ncols();
   
    double* A_ptr = A.raw();
    double* L_ptr = L.raw();
    for(int i=0;i<n;i++){
        int row = i*cols;
        double* A_row_i = A_ptr + row;
        double* L_row_i = L_ptr + row;

        double sum=0.0;
        for(int p=0;p<i;p++){
            sum += L_row_i[p]*L_row_i[p];
        }
        L_row_i[i] = std::sqrt(A_row_i[i] - sum);
        for(int j=0;j<i;j++){
            double sum=0.0;
            double* L_row_j = L_ptr + j*cols;
            for(int p=0;p<j;p++){
                sum += L_row_i[p] *L_row_j[p];
            }  
            L_row_i[j] = (A_row_i[j] - sum)/L_row_j[j];
        }
    }
    
}

inline void lu_solve(const Matrix& A, const std::vector<double>& b, std::vector<double>& x, std::string choice, bool print_parts = false)
{
    const int n = A.nrows();
    Matrix L(n, n);
    Matrix U(n, n);
    std::vector<int> P(n);
    x.resize(n);
    Matrix A_work = A;
    if(choice == "doolittle") lu_doolittle(A_work,L,U,P);
    else if(choice == "crout") lu_crout(A_work,L,U,P);
    else if(choice=="cholesky"){
        lu_cholesky(A_work,L);
        transpose(L,U);
    }
    else throw std::runtime_error("LU Options (doolittle/crout/cholesky)");
    if(print_parts){
        std::cout << "L Matrix:\n" << L << "\nU Matrix:\n" << U << std::endl;
    }

    std::vector<double> Pb(n);
    for (int i=0;i<n;i++) {
        Pb[i] = b[P[i]];
    }

    std::vector<double> y(n);
    forward_sub(L, Pb.data(), y.data());
    back_sub(U, y.data(), x.data());

}
// O(2 n^3)
inline void qr_gram_schmidt(
    const Matrix& A,
    Matrix& Q,
    Matrix& R
){
    const int rows = A.nrows();
    const int cols = A.ncols();
    const double* A_ptr = A.raw();
    double* Q_ptr = Q.raw();
    double* R_ptr = R.raw();

    Matrix V(rows, cols);
    double* V_ptr = V.raw();

    for(int k=0;k<cols;k++){
        for(int i=0;i<rows;i++){
            V_ptr[i*cols+k] = A_ptr[i*cols+k];
        }

        for(int j=0;j<k;j++){
            double dot_num = 0.0;
            double dot_den = 0.0;
            for(int i=0;i<rows;i++){
                dot_num += V_ptr[i*cols+j] * A_ptr[i*cols+k];
                dot_den += V_ptr[i*cols+j] * V_ptr[i*cols+j];
            }
            double scalar = dot_num / dot_den;
            for(int i=0;i<rows;i++){
                V_ptr[i*cols+k] -= scalar * V_ptr[i*cols+j];
            }
        }

        double v_norm = 0.0;
        for(int i=0;i<rows;i++){
            v_norm += V_ptr[i*cols+k] * V_ptr[i*cols+k];
        }
        v_norm = std::sqrt(v_norm);
        for(int i=0;i<rows;i++){
            Q_ptr[i*cols+k] = V_ptr[i*cols+k] / v_norm;
        }
    }

    for(int i=0;i<cols;i++){
        double* R_row_i = R_ptr + i*cols;
        for(int j=0;j<cols;j++){
            double sum = 0.0;
            for(int k=0;k<rows;k++){
                sum += Q_ptr[k*cols+i] * A_ptr[k*cols+j];
            }
            R_row_i[j] = sum;
        }
    }
}

// O(4/3 * n^3) each for Q and R
inline void qr_householder_inplace(
    Matrix& A,
    Matrix& QT
){
    const int n = A.nrows();
    const int cols = A.ncols();
    double* A_ptr = A.raw();
    double* QT_ptr = QT.raw();

    for(int i=0;i<n;i++){
        double* QT_ptr_i = QT_ptr + i*n;
        for(int j=0;j<n;j++){
            QT_ptr_i[j] = (i==j) ? 1.0 : 0.0;
        }
    }

    std::vector<double> u(n);
    for(int i=0;i<n-1;i++){
        double norm_Ai = 0.0;
        double* A_row_i = A_ptr + i*cols;
        for(int r=i;r<n;r++){
            norm_Ai += A_ptr[r*cols+i] * A_ptr[r*cols+i];
        }
        double Ai_cap = -(A_row_i[i] >= 0 ? 1.0 : -1.0) * std::sqrt(norm_Ai);
        
        for(int r=i;r<n;r++) u[r-i] = A_ptr[r*cols+i];
        u[0] = std::sqrt((Ai_cap - A_row_i[i]) / (2.0 * Ai_cap));
        for(int j=1;j<n-i;j++){
            u[j] = A_ptr[(i+j)*cols+i] / (-2.0 * Ai_cap * u[0]);
        }
        
        for(int c=i;c<n;c++){
            double dot = 0.0;
            for(int k=0;k<n-i;k++) dot += u[k] * A_ptr[(i+k)*cols+c];
            for(int r=0;r<n-i;r++) A_ptr[(i+r)*cols+c] -= 2.0 * u[r] * dot;
        }
        
        for(int c=0;c<n;c++){
            double dot = 0.0;
            for(int k=0;k<n-i;k++) dot += u[k] * QT_ptr[(i+k)*cols+c];
            for(int r=0;r<n-i;r++) QT_ptr[(i+r)*cols+c] -= 2.0 * u[r] * dot;
        }
    }
}

// O(2n^3) for each Q and R
inline void qr_givens_inplace(
    Matrix& R,
    Matrix& QT
){
    const int n = R.nrows();
    const int cols = R.ncols();
    double* R_ptr = R.raw();
    double* QT_ptr = QT.raw();

    for(int i=0;i<n;i++){
        double* QT_ptr_i = QT_ptr + i*n;
        for(int j=0;j<n;j++){
            QT_ptr_i[j] = (i==j) ? 1.0 : 0.0;
        }
    }

    for(int i=0;i<n-1;i++){
        for(int j=n-1;j>i;j--){
            int i1 = j;
            int i2 = j-1;
            double* R_row_i1 = R_ptr + i1*cols;
            double* R_row_i2 = R_ptr + i2*cols;
            double* QT_col_i1 = QT_ptr + i1*n;
            double* QT_col_i2 = QT_ptr + i2*n;
            if(std::abs(R_row_i1[i]) < 1e-12) continue;
            
            double a = R_row_i2[i];
            double b = R_row_i1[i];
            double r = std::sqrt(a*a + b*b);
            double s = b/r;
            double c = a/r;
            
            for(int k=0;k<cols;k++){
                double t1 = R_row_i2[k];
                double t2 = R_row_i1[k];
                R_row_i2[k] = c*t1 + s*t2;
                R_row_i1[k] = -s*t1 + c*t2;
            }
            for(int k=0;k<n;k++){
                double t1 = QT_col_i2[k];
                double t2 = QT_col_i1[k];
                QT_col_i2[k] = c*t1 + s*t2;
                QT_col_i1[k] = -s*t1 + c*t2;
            }
        }
    }
}

inline void qr_solve(
    const Matrix& A,
    std::vector<double>& b,
    std::vector<double>& x,
    std::string method = "householder"
){
    const int n = A.nrows();
    Matrix R = A;
    Matrix QT(n, n);
    x.resize(n);

    if(method == "householder") qr_householder_inplace(R, QT);
    else if(method == "givens") qr_givens_inplace(R, QT);
    else throw std::runtime_error("QR Options (householder/givens)");

    std::vector<double> y(n, 0.0);
    double* QT_ptr = QT.raw();
    for(int i=0;i<n;i++){
        y[i] = dot(QT_ptr + i*n, b.data(), n);
    }
    back_sub(R, y.data(), x.data());
}

// O(2n^2) per iteration
inline void jacobi_inplace(
    const Matrix& A,
    const double* _RESTRICT b,
    double* _RESTRICT x,
    double tol = 1e-6,
    int max_iter = 100
){
    const int n = A.nrows();
    const int cols = A.ncols();
    const double* A_ptr = A.raw();
    
    std::vector<double> x_k(n, 1.0);
    std::vector<double> x_k2(n, 0.0);
    int cycles = 0;

    for(int iter=0;iter<max_iter;iter++){
        cycles++;
        double max_diff = 0.0;
        for(int i=0;i<n;i++){
            double s = 0.0;
            const double* A_row_i = A_ptr + i*cols;
            for(int k=0;k<n;k++){
                if(k!=i) s += A_row_i[k] * x_k[k];
            }
            x_k2[i] = (b[i] - s) / A_row_i[i];
            double diff = std::abs(x_k2[i] - x_k[i]);
            if(diff > max_diff) max_diff = diff;
        }
        if(max_diff < tol){
            std::cout << "\nCycles : " << cycles << std::endl;
            for(int i=0;i<n;i++) x[i] = x_k2[i];
            return;
        }
        for(int i=0;i<n;i++) x_k[i] = x_k2[i];
    }
    std::cout << "\nCycles : " << cycles << std::endl;
    for(int i=0;i<n;i++) x[i] = x_k2[i];
}

inline void jacobi(const Matrix& A, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-6, int max_iter = 100){
    x.resize(A.nrows());
    jacobi_inplace(A, b.data(), x.data(), tol, max_iter);
}

// O(2n^2) per iter
inline void gauss_seidel_inplace(
    const Matrix& A,
    const double* _RESTRICT b,
    double* _RESTRICT x,
    double tol = 1e-6,
    int max_iter = 1000
){
    const int n = A.nrows();
    const int cols = A.ncols();
    const double* A_ptr = A.raw();
    
    for(int i=0;i<n;i++) x[i] = 1.0;
    std::vector<double> x_old(n);
    int cycles = 0;

    for(int iter=0;iter<max_iter;iter++){
        cycles++;
        for(int i=0;i<n;i++) x_old[i] = x[i];
        
        double max_diff = 0.0;
        for(int i=0;i<n;i++){
            double s = 0.0;
            const double* A_row_i = A_ptr + i*cols;
            for(int k=0;k<n;k++){
                if(k!=i) s += A_row_i[k] * x[k];
            }
            x[i] = (b[i] - s) / A_row_i[i];
            double diff = std::abs(x_old[i] - x[i]);
            if(diff > max_diff) max_diff = diff;
        }
        if(max_diff < tol){
            std::cout << "\nCycles : " << cycles << std::endl;
            return;
        }
    }
    std::cout << "\nCycles : " << cycles << std::endl;
}

inline void gauss_seidel(const Matrix& A, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-6, int max_iter = 1000){
    x.resize(A.nrows());
    gauss_seidel_inplace(A, b.data(), x.data(), tol, max_iter);
}

inline void successive_relaxation_inplace(
    double w,
    const Matrix& A,
    const double* _RESTRICT b,
    double* _RESTRICT x,
    double tol = 1e-6,
    int max_iter = 1000
){
    const int n = A.nrows();
    const int cols = A.ncols();
    const double* A_ptr = A.raw();
    
    for(int i=0;i<n;i++) x[i] = 1.0;
    std::vector<double> x_old(n);
    int cycles = 0;

    for(int iter=0;iter<max_iter;iter++){
        cycles++;
        for(int i=0;i<n;i++) x_old[i] = x[i];
        
        double max_diff = 0.0;
        for(int i=0;i<n;i++){
            double s = 0.0;
            const double* A_row_i = A_ptr + i*cols;
            for(int k=0;k<n;k++){
                if(k!=i) s += A_row_i[k] * x[k];
            }
            x[i] = (w*(b[i] - s))/A_row_i[i] + (1.0-w)*x_old[i];
            double diff = std::abs(x_old[i] - x[i]);
            if(diff > max_diff) max_diff = diff;
        }
        if(max_diff < tol){
            std::cout << "\nCycles : " << cycles << std::endl;
            return;
        }
    }
    std::cout << "\nCycles : " << cycles << std::endl;
}

inline void successive_relaxation(double w, const Matrix& A, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-6, int max_iter = 1000){
    x.resize(A.nrows());
    successive_relaxation_inplace(w, A, b.data(), x.data(), tol, max_iter);
}

// O(4n^2) per iter
inline void steepest_descent_inplace(
    const Matrix& A,
    const double* _RESTRICT b,
    double* _RESTRICT x,
    double tol = 1e-6,
    int max_iter = 100
){
    const int n = A.nrows();
    for(int i=0;i<n;i++) x[i] = 1.0;
    
    std::vector<double> r(n);
    std::vector<double> Ax(n);
    matvec(A, x, Ax.data());
    for(int i=0;i<n;i++) r[i] = b[i] - Ax[i];
    
    std::vector<double> Ar(n);

    for(int i=0;i<max_iter;i++){
        double r_dot = dot(r.data(), r.data(), n);
        if(std::sqrt(r_dot) < tol){
            std::cout << "Cycles : " << i << std::endl;
            return;
        }
        
        matvec(A, r.data(), Ar.data());
        double r_Ar = dot(r.data(), Ar.data(), n);
        double alpha = r_dot / r_Ar;
        scalar_mult_add(alpha,r.data(),x,n);
        scalar_mult_add(-alpha,Ar.data(),r.data(),n);

    }
    std::cout << "Cycles : " << max_iter << std::endl;
}

inline void steepest_descent(const Matrix& A, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-6, int max_iter = 100){
    x.resize(A.nrows());
    steepest_descent_inplace(A, b.data(), x.data(), tol, max_iter);
}

// O(2n^2) per iter
inline void conjugate_gradient_inplace(
    const Matrix& A,
    double* _RESTRICT b,
    double* _RESTRICT x,
    double tol = 1e-6,
    int max_iter = 100
){
    const int n = A.nrows();
    for(int i=0;i<n;i++) x[i] = 1.0;
    
    std::vector<double> r(n);
    std::vector<double> Ax(n);
    matvec(A, x, Ax.data());
    for(int i=0;i<n;i++) r[i] = b[i] - Ax[i];
    
    std::vector<double> d = r;
    std::vector<double> Ad(n);
    double b_norm = std::sqrt(dot(b, b, n));

    for(int i=0;i<max_iter;i++){
        double r_dot = dot(r.data(), r.data(), n);
        if(std::sqrt(r_dot)/(b_norm+1e-16) < tol){
            std::cout << "Cycles : " << i << std::endl;
            return;
        }
        
        double r_d = dot(r.data(), d.data(), n);
        matvec(A, d.data(), Ad.data());
        double d_Ad = dot(d.data(), Ad.data(), n);
        
        double alpha = r_d / d_Ad;
        scalar_mult_add(alpha,d.data(),x,n);
        scalar_mult_add(-alpha,Ad.data(),r.data(),n);

        
        double r_Ad = dot(r.data(), Ad.data(), n);
        double beta = -(r_Ad / d_Ad);
        if(std::abs(d_Ad) < 1e-14){
            std::cout << "near-zero denominator in CG beta calc\n";
            return;
        }
        for(int j=0;j<n;j++){
            d[j] = r[j] + beta * d[j];
        }
    }
    std::cout << "Cycles : " << max_iter << std::endl;
}

inline void conjugate_gradient(const Matrix& A,std::vector<double>& b, std::vector<double>& x, double tol = 1e-6, int max_iter = 100){
    x.resize(A.nrows());
    conjugate_gradient_inplace(A, b.data(), x.data(), tol, max_iter);
}