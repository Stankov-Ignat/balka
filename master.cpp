#include "head.h"
using namespace std;

// <============================================MAIN==============================================================>

int main (void)
{
  double data[16];    // [x_k, x_q1, x_q2, x_r1, x_r2, x_r3, x_m, x_z, x_p][q1, q2, p, k, m, E, J]

  read_data ("data.txt", data);

//--------------------------------------------------------------------------------- создание сетки   
/* ебло ты утиное */

  int L = 25;             // Матрица A будет размера 2L x 2L.      L - 1 = колво конеч эл.
  
  double* x_mesh = new double[L];

  int index_v[4];      // x_r1  x_r2  x_r3  x_z / Индексы узлов в глобальном смысле

  if (make_mesh (L, data, x_mesh, index_v))    return -1; 
  
  //double* x_mesh = make_uniform_mesh (L, index_v, 0);

//---------------------------------------------------------------------------------- выделение памяти под A и F 

  double** A = new double*[2 * L];                                                   
  for (size_t i = 0; i < 2 * L; i++)   
    { 
      A[i] = new double[2 * L];     fill (A[i], A[i] + 2 * L, 0.0);
    } 

  double* F = new double[2 * L];    fill (F, F + 2 * L, 0.0);  

//-----------------------------------------------------------------------------------

  for (int k = 0; k < L - 1; k++)               // k - индекс конечного элемента
    {
      for (int i = 0; i < 4; i++)
        {
          for (int j = 0; j < 4; j++)
            {
              A[2 * k + i][2 * k + j] += a_res (p_def (i+1, j+1), x_mesh[k], x_mesh[k + 1], data[14], data[15], data[0], data[12]);
            }

          F[2 * k + i] += f_res (i + 1, x_mesh[k], x_mesh[k + 1], data[9], data[10], data[1], data[2], data[13], data[6], data[11], data[8]);
        }
    }

//-------------------------------------------------------------------------------------

  double *solve_conju_ = new double[2 * L];    fill (solve_conju_, solve_conju_ + 2 * L, 0.0);

  double react[5];    // r1 r2 r3 r_z m_z   |  последнее число + 2*m = m_z    опытным путем выяснил
  
  conjugate_gradient_s_limit (2 * L, A, F, solve_conju_, index_v, react);
  
  //cout << "x_mesh = ";    show_vect (x_mesh, L, 1);
  cout << "reactions = ";     show_vect (react, 5, 1);
  //cout << "solve_limit = ";   show_vect (solve_conju_, 2 * L, 1);

  array_to_file("solve.txt", solve_conju_, 2 * L);

//-------------------------------------------------------------------------------------
/*
  double** A_tmp = new double*[2 * L];                                                   
  for (size_t i = 0; i < 2 * L; i++) { A_tmp[i] = new double[2 * L]; } 

  double* F_tmp = new double[2 * L];

  for (size_t i = 0; i < 2 * L; i++)
    {
      for (size_t j = 0; j < 2 * L; j++)
        {
          A_tmp[i][j] = A[i][j];
        }
      F_tmp[i] = F[i];
    }

  int size_A = 2 * L;
  
  rm_index_v (A_tmp, F_tmp, size_A, index_v);         // удаление части строк и столбцов

  size_t del_count_A = size_A;

  double *solve_conju = new double[2 * L];    fill (solve_conju, solve_conju + 2 * L, 0.0);
  
  conjugate_gradient (size_A, A_tmp, F_tmp, solve_conju);
  
  verify_solution (A_tmp, solve_conju, F_tmp, size_A, 1.e-13, false);

  add_v_to_solve (solve_conju, size_A, index_v);

  auto l2_norm = [](double* a, double*b, int size) {double res = 0.0;   for (size_t i = 0; i < size; i++) res +=  pow (a[i] - b[i], 2); return sqrt (res); };

  cout << "Conju_grad_2_ways: l2_norm = " << l2_norm (solve_conju, solve_conju_, 2 * L) << endl;
*/
//-----------------------------------------------СЕТКА ДЛЯ ГРАФИКИ----------------------------------------

  int show_count = 1e3;
  double* show_mesh = new double[show_count];  // сетка на [0, 25] для отображдения из локальных координат

  double* show_w = new double[show_count];         // прогиб
  double* show_theta = new double[show_count];     // угол производная прогиба
  double* show_moment = new double[show_count];    // угол производная прогиба
  double* show_q = new double[show_count];         // угол производная прогиба

  make_solution (show_count, show_mesh, x_mesh, solve_conju_, show_w, show_theta, show_moment, show_q, L, data);

  write_arrays_to_file("output.txt", show_mesh, show_w, show_theta, show_moment, show_q, show_count);

  system("sage main.sage");

//---------------------------------------------------------------------------------------

  delete[] x_mesh;
  delete[] show_mesh;
  delete[] show_w;
  delete[] show_theta;
  delete[] show_moment;
  delete[] show_q;

  delete[] F;
 
  for (size_t i = 0; i < 2 * L; i++)
    {
      if (A[i] != nullptr) 
        {
          delete[] A[i];
          A[i] = nullptr; // Обнуление указателя
        }
    }
  delete[] A;
  
  return 0;
}