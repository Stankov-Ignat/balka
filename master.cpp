#include "head.h"
using namespace std;

/* Функция для решения системы линейных уравнений методом сопряженных градиентов */
void conjugate_gradient_s_limit (int size, double **A, double *B, double *x_conju, int* index_v, double* react,
                                 int iter_max = 1e2, double epsilon = 1e-14, bool showing_param = 0)
{
  double *p = new double[size];    // p_k направление спуска 
  double *r = new double[size];    // r_k невязка
  double *tt = new double[size];   /* tt = A * p_k */
  double alpha, beta, tmp1, nevyzka, B_norm = 0;
  
  // возвращает 1 если позиция i в блок - векторе 2*L в тех индексах, которые надо занулить
  auto is_in_index_v = [](int i, int* index_v) 
     { return (i == 2 * index_v[0] || i == 2 * index_v[1] || i == 2 * index_v[2] || i == 2 * index_v[3] || i == (2 * index_v[3] + 1));}; 

  auto cond_vect = [](double* a, int size, int* index_v) {
      for (size_t i = 0; i < 4; i++)  {a[2 * index_v[i]] = 0.0;}
      a[2 * index_v[3] + 1] = 0.0;
    }; 

  cond_vect (x_conju, size, index_v);

  /* Инициализация r_0 и p_0 */
  for (int i = 0; i < size; i++)
    {
      r[i] = B[i];
      for (size_t k = 0; k < size; k++)  { r[i] -= A[i][k] * x_conju[k]; }  // случай ненулевого начального условия x_0
      p[i] = r[i];
      B_norm += B[i] * B[i];  /* Норма B */
    }

  cond_vect (p, size, index_v);
  B_norm = sqrt(B_norm);

  /* Итерации метода сопряженных градиентов */
  for (int count = 0; count < iter_max; count++)  /* Счетчик итераций */
    {
      alpha = 0.0;
      tmp1 = 0.0;
      beta = 0.0;

      /* Вычисление A * p_k */
      for (int i = 0; i < size; i++)
        {
          tt[i] = 0.0;
          if (!is_in_index_v (i, index_v))
            {
              for (int k = 0; k < size; k++)  /* Вычисляем tt = A * p_k */
                {
                  tt[i] += A[i][k] * p[k];
                }
            }
          tmp1 += r[i] * r[i];    /* tmp1 = (r_k, r_k) */
          alpha += p[i] * tt[i];  /* alpha = (p_k, A * p_k) */
        }

            /* Проверка деления на ноль */
      if (fabs(alpha) < 1e-12)
        {
          cout << "Ошибка: деление на ноль в alpha" << endl;
          break;
        }

      /* Вычисление alpha */
      alpha = tmp1 / alpha;

      /* Обновление x_k+1 и r_k+1 */
      for (int i = 0; i < size; i++)
        {
          x_conju[i] += alpha * p[i];
          r[i] -= alpha * tt[i];
        }

      /* Вычисление beta */
      for (int i = 0; i < size; i++)
        {
          beta += r[i] * r[i];
        }

      /* Норма остатка */
      nevyzka = beta;
      beta = beta / tmp1;

      /* Метод Флетчера - Ривса нужен рестарт на каждом size+1 иттерации нужно забывать направление спуска */
      if (count % (size + 1) == 0)
        {
          beta = 0.0;
        }

      /* Обновление p_k+1 */
      for (int i = 0; i < size; i++)
        {
          p[i] = r[i] + beta * p[i];
        }
      cond_vect (p, size, index_v);

      /* Проверка условия останова */
      if (nevyzka < epsilon * B_norm)
        {
          if (showing_param)
            {
              cout << "Метод сопряженных градиентов завершился после " << count + 1 << " итераций." << endl;
            }
          break;
        }
    }
  // after itterations
  for (int i = 0; i < size; i++)
    {
      r[i] = -B[i];
      for (size_t k = 0; k < size; k++)  { r[i] += A[i][k] * x_conju[k]; }  // случай ненулевого начального условия x_0
    }
  //cout << "r = ";   show_vect (r, size, 0);
  
  for (size_t i = 0; i < 4; i++)
    {
      react[i] = r[2 * index_v[i]];
    }
  react[4] = r[2 * index_v[3] + 1];

  delete[] p;
  delete[] r;
  delete[] tt;
}

// Функция для поиска совпадений
void find_indices(double* x_mesh, int mesh_size, double* values, int* found_indices, double tolerance = 1e-6) 
{
  int values_size = 4;

    // Проход по массиву значений
  for (int i = 0; i < values_size; ++i) 
    {
      for (int j = 0; j < mesh_size; ++j) 
        {
          if (abs(x_mesh[j] - values[i]) < tolerance) 
            { 
              found_indices[i] = j;
              break; // Найдено, переходим к следующему значению
          }
      }
  }
}

// <==============================================================================================================>








// <============================================MAIN==============================================================>
int main (void)
{
  double data[16];    // [x_k, x_q1, x_q2, x_r1, x_r2, x_r3, x_m, x_z, x_p][q1, q2, p, k, m, E, J]
  read_data ("data.txt", data);

//--------------------------------------------------------------------------------- создание сетки   
/* ебло ты утиное */
/*
  int L = 20;             // Матрица A будет размера 2L x 2L.      L - 1 = колво конеч эл.
  double* x_mesh = new double[L];

  int index_v[4];      // x_r1  x_r2  x_r3  x_z / Индексы узлов в глобальном смысле
  index_v[3] = L - 1;     

  if (make_mesh(L, data, x_mesh, index_v))
    return -1; 

  array_to_file("mesh.txt", x_mesh, L);
*/
//---------------------------------------------------------------------------------- ручной режим управления сосать америка

  int n = 0;
  int L = 25 * static_cast<int>(pow(2, n)) + 1;

  double h = 1 / pow (2, n);
  double* x_mesh = new double[L];
  for (int i = 0; i < L; i++)
    x_mesh[i] = i * h;
  
  int index_v[4];

  double values[4] = {4, 12, 16, 25};

  find_indices(x_mesh, L, values, index_v); 

  show_vect_int (index_v, 4, 1);
  
  array_to_file("mesh.txt", x_mesh, L);

//---------------------------------------------------------------------------------- выделение памяти под A и F 

  double** A = new double*[2 * L];                                                   
  for (size_t i = 0; i < 2 * L; i++)   
    { 
      A[i] = new double[2 * L];     fill (A[i], A[i] + 2 * L, 0.0);
    } 

  double* F = new double[2 * L];    fill (F, F + 2 * L, 0.0);  

  double** A_tmp = new double*[2 * L];                                                   
  for (size_t i = 0; i < 2 * L; i++) { A_tmp[i] = new double[2 * L]; } 

  double* F_tmp = new double[2 * L];      

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

  //cout << "\nA и F:\n"; 
  //show_array (A, 2 * L);
  //show_vect (F, 2 * L, 1);

  for (size_t i = 0; i < 2 * L; i++)
    {
      for (size_t j = 0; j < 2 * L; j++)
        {
          A_tmp[i][j] = A[i][j];
        }
      F_tmp[i] = F[i];
    }

//-------------------------------------------------------------------------------------

  double *solve_conju_ = new double[2 * L];    fill (solve_conju_, solve_conju_ + 2 * L, 0.0);

  double react[5];
  
  conjugate_gradient_s_limit (2 * L, A, F, solve_conju_, index_v, react);
  
  cout << "x_mesh = ";    show_vect (x_mesh, L, 1);
  cout << "react = ";     show_vect (react, 5, 1);
  //cout << "solve_1 = ";   show_vect (solve_conju_, 2 * L, 1);

  verify_solution (A, solve_conju_, F, 2 * L);
  
  array_to_file("solve.txt", solve_conju_, 2 * L);

//-------------------------------------------------------------------------------------
/*
  int size_A = 2 * L;

  rm_index_v (A_tmp, F_tmp, size_A, index_v);         // удаление части строк и столбцов

  size_t del_count_A = size_A;

//----------------------------------------------------------------------------------------

  double *solve_conju = new double[2 * L];    fill (solve_conju, solve_conju + 2 * L, 0.0);

  conjugate_gradient (size_A, A_tmp, F_tmp, solve_conju);
  verify_solution (A_tmp, solve_conju, F_tmp, size_A);

  add_v_to_solve (solve_conju, size_A, index_v);
  //cout << "solve_2 = ";   show_vect (solve_conju, 2 * L, 1);
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
 /*
  for (size_t i = 0; i < 2 * L; i++)
    {
      if (A[i] != nullptr) 
        {
          delete[] A[i];
          A[i] = nullptr; // Обнуление указателя
        }
    }
  delete[] A;
  */
  return 0;
}