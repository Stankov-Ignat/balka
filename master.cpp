#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <stdexcept> // Для runtime_error

using namespace std;

//---------------------------------------------------------------------------------------------------------

void show_vect (double*, int, int);

void show_array (double**, int);

/* Функция для решения системы линейных уравнений методом сопряженных градиентов */
void conjugate_gradient (int size, double **A, double *B, double *x_conju, int iter_max = 1e2, double epsilon = 1e-14, bool showing_param = 0)
{
  double *p = new double[size];    /* p_k */
  double *r = new double[size];    /* r_k */
  double *tt = new double[size];   /* tt = A * p_k */
  double alpha, beta, tmp1, nevyzka, B_norm = 0;

  /* Инициализация r_0 и p_0 */
  for (int i = 0; i < size; i++)
    {
      r[i] = B[i];
      for (size_t k = 0; k < size; k++)  { r[i] -= A[i][k] * x_conju[k]; }  // случай ненулевого начального условия x_0
      p[i] = r[i];
      B_norm += B[i] * B[i];  /* Норма B */
    }

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
          for (int k = 0; k < size; k++)  /* Вычисляем tt = A * p_k */
            {
              tt[i] += A[i][k] * p[k];
            }

          tmp1 += r[i] * r[i];  /* tmp1 = (r_k, r_k) */
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

      // cout << "Итерация " << count << " Невязка: " << sqrt(beta) << endl;

      /* Норма остатка */
      nevyzka = beta;
      beta = beta / tmp1;

      /* Еще он звисит но начальной точки */
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

  delete[] p;
  delete[] r;
  delete[] tt;
}

void LU_(int size, double **A, double *B, double* x_LU)
{
  double** L = new double *[size];                                                                               
  double** U = new double *[size];        

  for (size_t i = 0; i < size; i++)
    {
      U[i] = new double[size];   L[i] = new double[size]; 
    }        

  for (size_t i = 0; i < size; i++) 
   {
      /* Вычисление элементов U для строки i */
      for (size_t j = i; j < size; j++) 
        {
          U[i][j] = A[i][j];
          for (size_t k = 0; k < i; k++) 
            {
              U[i][j] -= L[i][k] * U[k][j];
            }
        }

      /* Вычисление элементов L для столбца i */
      for (size_t j = i + 1; j < size; j++) 
        {
          L[j][i] = A[j][i];
          for (size_t k = 0; k < i; k++) 
            {
              L[j][i] -= L[j][k] * U[k][i];
            }
          L[j][i] /= U[i][i];
        }

      /* Диагональный элемент L[i][i] всегда равен 1 */
      L[i][i] = 1.0;
  }
  
  for (int i = 0; i < size; i++)                   // обратный ход L
    {
      x_LU[i] = B[i];

      for (size_t j = 0; j < i; j++)
        {
          x_LU[i] -= L[i][j] * x_LU[j];
        }
    }

  for (int i = size - 1; i > -1; i--)             // обратный ход U
    {
      for (size_t j = i + 1; j < size; j++)
        {
          x_LU[i] -= U[i][j] * x_LU[j];
        }
      x_LU[i] /= U[i][i];
    }

  for (size_t i = 0; i < size; i++) 
    {
      delete[] L[i];
      delete[] U[i];
    }
  delete[] L;   delete[] U;  
}


/* базисные функции в порядке {N1, L1, N2, L2} */
double base_func (int p, double t)
{
  switch (p)
    {
      case 1: return 1 - 3 * pow(t, 2) + 2 * pow(t, 3);           // N1                           
      case 2: return t - 2 * t * t + pow(t, 3);                   // L1
      case 3: return 3 * pow(t, 2) - 2 * pow(t, 3);               // N2
      case 4: return -t * t + pow(t, 3);                          // L2
      default: 
        {
          cout << "OUT OF RANGE base_func" << endl;   return 0.0;
        }
    }
  return 0;
}

/* производные базисных функций */
double base_d_func (int p, double t)
{
  switch (p)
    {
      case 1: return -6 * t + 6 * pow(t, 2);                  // dN1                           
      case 2: return 1 - 4 * t + 3 * pow(t, 2);               // dL1
      case 3: return 6 * t - 6 * pow(t, 2);                   // dN2
      case 4: return -2 * t + 3 * pow(t, 2);                  // dL2
      default: 
        {
          cout << "OUT OF RANGE base_func" << endl;   return 0.0;
        }
    }
  return 0;
}

double base_dd_func (int p, double t)
{
  switch (p)
    {
      case 1: return -6 + 12 * t;              // ddN1                           
      case 2: return -4 + 6 * t;               // ddL1
      case 3: return  6 - 12 * t;              // ddN2
      case 4: return -2 + 6 * t;               // ddL2
      default: 
        {
          cout << "OUT OF RANGE base_func" << endl;   return 0.0;
        }
    }
  return 0;
}

double base_ddd_func (int p, double t)
{
  switch (p)
    {
      case 1: return  12;              // dddN1                           
      case 2: return   6;              // dddL1
      case 3: return -12;              // dddN2
      case 4: return   6;              // dddL2
      default: 
        {
          cout << "OUT OF RANGE base_func" << endl;   return 0.0;
        }
    }
  return 0;
}

/* меняет i, j на номера базисных функций по которым считался интеграл в a_base */
void p_change (int p, int &i, int &j)
{
  switch (p)
    {
      case 0: i = 1; j = 1;   break; 
      case 1: i = 1; j = 2;   break;
      case 2: i = 1; j = 3;   break;
      case 3: i = 1; j = 4;   break;
      case 4: i = 2; j = 2;   break;
      case 5: i = 2; j = 3;   break;
      case 6: i = 2; j = 4;   break;
      case 7: i = 3; j = 3;   break;
      case 8: i = 3; j = 4;   break;
      case 9: i = 4; j = 4;   break;
      default: 
        {
          cout << "OUT OF RANGE p_change" << endl;  break;
        }
    }
}

/* по i,j --> p номер в a_base */
int p_def(int i, int j) 
{
  int p = -1;
  if (i > j)  swap(i, j); // i <= j

  if (i == 1 && j == 1) return 0;  // 11
  if (i == 1 && j == 2) return 1;  // 12
  if (i == 1 && j == 3) return 2;  // 13
  if (i == 1 && j == 4) return 3;  // 14
  if (i == 2 && j == 2) return 4;  // 22
  if (i == 2 && j == 3) return 5;  // 23
  if (i == 2 && j == 4) return 6;  // 24
  if (i == 3 && j == 3) return 7;  // 33
  if (i == 3 && j == 4) return 8;  // 34
  if (i == 4 && j == 4) return 9;  // 44  
  
  return p; 
}

/* на каждои отрезке неопределенный инт a_base от 0 до 1 */
double a_base (int p)
{
  switch (p)
    {
      case 0: return  12;      // 36 * t - 72 * pow(t, 2) + 48 * pow(t, 3);             // 11
      case 1: return   6;      // 6 * (4 * t - 7 * pow(t, 2) + 4 * pow(t, 3));          // 12
      case 2: return -12;      // -36 * t + 72 * pow(t, 2) - 48 * pow(t, 3);            // 13
      case 3: return   6;      // 6 * (2 * t - 5 * pow(t, 2) + 4 * pow(t, 3));          // 14
      case 4: return   4;      // 16 * t - 24 * pow(t, 2) + 12 * pow(t, 3);             // 22
      case 5: return  -6;      // -6 * (4 * t - 7 * pow(t, 2) + 4 * pow(t, 3));         // 23
      case 6: return   2;      // 2 * (4 * t - 9 * pow(t, 2) + 6 * pow(t, 3));          // 24
      case 7: return  12;      // 36 * t - 72 * pow(t, 2) + 48 * pow(t, 3);             // 33
      case 8: return  -6;      // -6 * (2 * t - 5 * pow(t, 2) + 4 * pow(t, 3));         // 34
      case 9: return   4;      // 4 * t - 12 * pow(t, 2) + 12 * pow(t, 3);              // 44
      default: 
        {
          cout << "OUT OF RANGE a_base" << endl;   return 0.0;
        }
    }
  return 0;
}

/* считает определенный интеграл от x_i до x_(i+1) и добавляет работу пружины если пружина на K */
double a_res (int p, double x_i, double x_right, double E, double J, double x_k, double k)
{
  double h_cur = x_right - x_i;
  double res = E * J * pow(1 / h_cur, 3) * a_base(p);  // 1/h^3 вылезает из вторых производных и якобиана преобразования

/* 
  cout << "h_cur = " << h_cur << "  ";
  cout << "p = " << p << "  x_i = " << x_i << "  x_i+1 = " << x_right << "  res = " << res << "   int(1) = " << a_base(p);
  cout << "   E*J = " << E*J << "    1/h^3 = " << pow(1 / h_cur, 3) << endl;
*/

  if (x_k < x_i || x_k > x_right)
    {
      return res;
    }
  else /* отображаем элемент, на который попал x_k в [0,1] получаем кси_к (определенную от 0 до 1) и на ней считаем произведение */
    {
      int i,j;    p_change(p, i, j);
      //cout << "  add = " << k * base_func(i, (x_k - x_i) / h_cur) * base_func(j, (x_k - x_i) / h_cur) << endl;
      return res - k * base_func(i, (x_k - x_i) / h_cur) * base_func(j, (x_k - x_i) / h_cur);
    }  
}

/* считает определенный интеграл он B_i * phi_j, где B_i линейный базис */
double f_integral (int i, int j)
{
  if (i == 1 && j == 1)   return 3 / 20;
  if (i == 1 && j == 2)   return 1 / 30;
  if (i == 1 && j == 3)   return 7 / 20;
  if (i == 1 && j == 4)   return -1 / 20;

  if (i == 2 && j == 1)   return 7 / 20;
  if (i == 2 && j == 2)   return 1 / 20;
  if (i == 2 && j == 3)   return 3 / 20;
  if (i == 2 && j == 4)   return -1 / 30;

  return 0.0;
}

double q_global (double x, double q1, double q2, double x_q1, double x_q2)
{
  if (x < x_q1 || x > x_q2)   return 0.0;

  return q1 + (q2 - q1) / (x_q2 - x_q1) * (x - x_q1);
}


/* подсчет <f,phi_p> с распределенной нагрузкой, моментом и сосредоточенной силой 
   p - индекс базисных эрмитовых от 1 до 4 */
double f_res (int p, double x_i, double x_right, double q1, double q2, double x_q1, double x_q2, double m, double x_m, double p_force, double x_p)
{
  double res = 0.0;

/* часть с распределенной нагрузкой q 
   по построению сетки конечный элемент может либо лежать целиком на распред. нагрузке, либо целиком не лежать.
   поэтому проверка лежит ли середина правее x_i */
  double median = (x_i + x_right) / 2;
  if (median > x_q1 && median < x_q2)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   q_raspred" << endl;
      double q1_, q2_;                              // значения q(x) на концах конечного элемента
      q1_ = q_global(x_i, q1, q2, x_q1, x_q2);
      q2_ = q_global(x_right, q1, q2, x_q1, x_q2);
      res += (x_right - x_i) * (q2_ * f_integral(1, p) + q1_ * f_integral(2, p));
    }


/* часть с моментом m */  
  if (x_m >= x_i && x_m <= x_right)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   moment" << endl;
      res += m * base_d_func(p, (x_m - x_i) / (x_right - x_i));
    }

/* часть с силой p. Переведенная в локальные координаты */  
  if (x_p >= x_i && x_p <= x_right)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   force" << endl;
      res -= p_force * base_func(p, (x_p - x_i) / (x_right - x_i));
    }
  return res;  
}

void show_array (double** array, int size)
{
  cout << endl;
  for (int i = 0; i < size; ++i) 
    {
      for (int j = 0; j < size; ++j) 
        {
          cout << setw(10) << setprecision(3) << array[i][j] << " ";
          //cout << setw(12) << array[i][j] << " ";
        }
      cout << endl;
    }
    cout << endl;
}

void show_vect (double* vect, int size, int mode)
{
  for (int i = 0; i < size; ++i) 
    {
      //cout << setw(8) << setprecision(3) << vect[i] << "\n";
      cout << setprecision(4) << vect[i] << "  ";
      if (mode == 0)
        {
          cout << endl;
        }
    }
  cout << endl;  
}

void sort(double* array, int size, int* index_v) 
{
  int index_count = 3;
  int original_indices[size];

  for (int i = 0; i < size; i++) 
    original_indices[i] = i;
   
  for (int i = 0; i < size; i++)
  {
    for (int j = i + 1; j < size; j++)
      {
        if (array[i] > array[j])
          {
            swap(array[i], array[j]);
            swap(original_indices[i], original_indices[j]);
          }
      }
  }

  // Обновление массива index_v на основе новых индексов
  for (int i = 0; i < index_count; i++) 
    {
      int original_index = static_cast<int>(index_v[i]);
      
      for (int j = 0; j < size; j++)
        {
          if (original_indices[j] == original_index) 
            {
              index_v[i] = j; // Новый индекс после сортировки
              break;
            }
        }
    }
}

/* 7 обязательных узлов: начало, конец, x_q1, x_r1, x_r2, x_q2, x_r3 */
int make_mesh (int L, double* data, double* x_mesh, int* index_v) 
{
  if (L < 7)
    {
      cout << "ERROR: L < 7\n";  return -1; 
    }

  double h = 25.0 / double(L - 1 - 5);      // равномерный шаг с выбросом обязательных узлов

  for (int i = 0; i < L - 5; i++)
    {
      x_mesh[i] = h * i;
    }

  // [x_k, x_q1, x_q2, x_r1, x_r2, x_r3, x_m, x_z, x_p][q1, q2, p, k, m, E, J]
  x_mesh[L - 5] = data[1]; /* x_q1 */
  x_mesh[L - 4] = data[2]; /* x_q2 */

  x_mesh[L - 3] = data[3]; /* x_r1 */
  x_mesh[L - 2] = data[4]; /* x_r2 */
  x_mesh[L - 1] = data[5]; /* x_r3 */

  for (size_t i = 0; i < 3; i++)  // x_r1  x_r2  x_r3
    {
      index_v[i] = L - 3 + i;
    }

  sort (x_mesh, L, index_v);
  return 0;
}

void read_data (const string& filename, double* data)
{
  ifstream inputfile(filename);
  string line;    double number = 0;
  for (int i = 0; i < 16; i++)
    {
      if (getline(inputfile, line))
        {
          istringstream iss(line);
          if (iss >> number)  data[i] = number;
        }
    }
  inputfile.close();
}

// Функция для проверки корректности решения
bool verify_solution(double **A, double *Q, double *F, int size, double tolerance = 1e-8, bool showing_param = 0)
{
  // Вычисляем A * Q
  double *computed_F = new double[size];
  for (int i = 0; i < size; i++)
    {
      computed_F[i] = 0.0;
      for (int j = 0; j < size; j++)
        {
          computed_F[i] += A[i][j] * Q[j];
        }
    }

  for (int i = 0; i < size; i++)  {computed_F[i] = F[i] - computed_F[i]; }
  //out << "B - A*X = ";   show_vect (computed_F, size, 1);

  // Вычисляем норму разности между F и A * Q
  double norm_diff = 0.0;
  for (int i = 0; i < size; i++)
    {
      norm_diff += computed_F[i] * computed_F[i];
    }
  norm_diff = sqrt(norm_diff);

  // Освобождаем память
  delete[] computed_F;

  // Проверяем, укладывается ли норма в допустимое отклонение
  if (norm_diff < tolerance)
    {
      if (showing_param)
        {
          cout << "Solution verification success: l2_abs norm diff = " << norm_diff << endl;
        }
      return true; // Решение корректно
    }
  else
    {
      if (showing_param)
        {
          cout << "Solution verification failed: l2_abs norm diff = " << norm_diff << endl;
        }
      return false; // Решение некорректно
    }
}

inline void rm_row_column (double** A, double* F, int& size, int rm_index) 
{
  // Удаляем столбец, сдвигая элементы внутри каждой строки
  for (size_t i = 0; i < size; i++) 
    {
      for (size_t j = rm_index; j < size - 1; j++) 
        {
          A[i][j] = A[i][j + 1]; // Сдвиг элементов влево
        }
    }

  for (size_t i = rm_index; i < size - 1; i++)    // Сдвиг указателей на строки
    {
      F[i] = F[i + 1];
      A[i] = A[i + 1]; 
    }

  //delete[] A[size - 1];

  size--;
}

/* index_v = [x_r1, x_r2, x_r3, x_z] */
void rm_index_v (double** A, double* F, int& size, int* index_v)
{
  rm_row_column (A, F, size, 2 * index_v[3] + 1);   // rm last row and column

  for (int i = 3; i > - 1; i--)
    {
      rm_row_column (A, F, size, 2 * index_v[i]);
    }
}

/* добавляет 0 в массив на место index || size = размер укороченного массива */
inline void add_num_to_array (double* array, int& size, int index)
{
  double tmp = array[index];
  array[index] = 0.0;

  for (size_t i = index + 1; i < size + 1; i++)
    {
      swap (tmp, array[i]);
    }

  size++;  
}

inline void add_v_to_solve (double* array, int& size, int* index_v)
{
  for (size_t i = 0; i < 4; i++)
    {
      add_num_to_array (array, size, 2 * index_v[i]);
    }
  add_num_to_array (array, size, 2 * index_v[3] + 1);  
}

void write_arrays_to_file(const char* filename, double* show_mesh, double* show_w, double* show_theta, 
                          double* show_moment, double* show_q, int show_count)
{
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cout << "Ошибка: не удалось открыть файл для записи!" << std::endl;
        return;
    }

    int width = 15;
    for (int i = 0; i < show_count; i++) {
        outfile << scientific << setw (width) << show_mesh[i] 
                << setw (width) << show_w[i] << " " 
                << setw (width) << show_theta[i] << " " 
                << setw (width) << show_moment[i] << " " 
                << setw (width) << show_q[i] << "\n";
    }

    outfile.close(); 
    //cout << "Данные успешно записаны в файл: " << filename << endl;
}

// <==============================================================================================================>








// <============================================MAIN==============================================================>
int main (void)
{
  double data[16];    // [x_k, x_q1, x_q2, x_r1, x_r2, x_r3, x_m, x_z, x_p][q1, q2, p, k, m, E, J]
  read_data ("data.txt", data);

//--------------------------------------------------------------------------------- создание сетки   
/* ебло ты утиное */

  int L = 7;             // Матрица A будет размера 2L x 2L.      L - 1 = колво конеч эл.
  double* x_mesh = new double[L];
  int index_v[4];      // x_r1  x_r2  x_r3  x_z   Индексы узлов в глобальном смысле
  index_v[3] = L - 1;     

  if (make_mesh(L, data, x_mesh, index_v))
    return -1; 

//---------------------------------------------------------------------------------- выделение памяти под A и F 

  double** A = new double*[2 * L];                                                   
  for (size_t i = 0; i < 2 * L; i++)   
    { 
      A[i] = new double[2 * L];     memset(A[i], 0.0, 2 * L * sizeof(double));
    } 

  double* F = new double[2 * L];    memset(F, 0.0, 2 * L * sizeof(double));  

  double** A_mini = new double*[4];                                                   
  for (size_t i = 0; i < 4; i++)   
    { 
      A_mini[i] = new double[4];
    } 

  double* F_mini = new double[4];   

//-----------------------------------------------------------------------------------

  for (int k = 0; k < L - 1; k++)               // k - индекс конечного элемента
    {
      for (int i = 0; i < 4; i++)
        {
          for (int j = 0; j < 4; j++)
            {
              A[2 * k + i][2 * k + j] += a_res(p_def(i+1, j+1), x_mesh[k], x_mesh[k + 1], data[14], data[15], data[0], data[12]);
              //A_mini[i][j] = a_res(p_def(i+1, j+1), x_mesh[k], x_mesh[k + 1], data[14], data[15], data[0], data[12]);
            }

          F[2 * k + i] += f_res (i + 1, x_mesh[k], x_mesh[k + 1], data[9], data[10], data[1], data[2], data[13], data[6], data[11], data[8]);
          //cout << "k = " << k + 1 << endl;
          //F_mini[i] = f_res (i + 1, x_mesh[k], x_mesh[k + 1], data[9], data[10], data[1], data[2], data[13], data[6], data[11], data[8]);
        }
      //show_array(A_mini, 4);
      //show_vect(F_mini, 4, 0); cout<<endl;
    }

  //cout << "\nA и F:\n"; 
  //show_array(A, 2 * L);
  //show_vect(F, 2 * L, 0);

//-------------------------------------------------------------------------------------

  int size_A = 2 * L;
  rm_index_v (A, F, size_A, index_v);         // удаление части строк и столбцов

//----------------------------------------------------------------------------------------

  double *solve_lu = new double[2 * L];       fill (solve_lu, solve_lu + 2 * L, 0.0);
  double *solve_conju = new double[2 * L];    fill (solve_conju, solve_conju + 2 * L, 0.0);

  LU_ (size_A, A, F, solve_lu);
  verify_solution (A, solve_lu, F, size_A);

  conjugate_gradient (size_A, A, F, solve_conju, 1e2);
  verify_solution (A, solve_conju, F, size_A);

  //add_v_to_solve (solve_lu, size_A, index_v);
  add_v_to_solve (solve_conju, size_A, index_v);
  
  //cout << "\nQ_conj: ";    show_vect (solve_conju, size_A, 1);

//-----------------------------------------------СЕТКА ДЛЯ ГРАФИКИ----------------------------------------

  int show_count = 1e2;
  double* show_mesh = new double[show_count];  // сетка на [0, 25] для отображдения из локальных координат
  
  double show_h = 25.0 / (show_count - 1);
  for (size_t i = 0; i < show_count; i++)
    {
      show_mesh[i] = i * show_h;
    }

  double* show_w = new double[show_count];         // прогиб
  double* show_theta = new double[show_count];     // угол производная прогиба
  double* show_moment = new double[show_count];    // угол производная прогиба
  double* show_q = new double[show_count];         // угол производная прогиба

  int count = 0;
  for (int i = 0; i < L - 1; i++)   // по конечным элементам  [x_mesh[i], x_mesh[i + 1]]
    {
      int h_current = x_mesh[i + 1] - x_mesh[i];

      while (show_mesh[count] >= x_mesh[i] && show_mesh[count] < x_mesh[i + 1])
        {
          show_w[count] = solve_conju[2 * i] * base_func (1, (show_mesh[count] - x_mesh[i]) / h_current) +
                          solve_conju[2 * i + 1] * base_func (2, (show_mesh[count] - x_mesh[i]) / h_current) +
                          solve_conju[2 * (i + 1)] * base_func (3, (show_mesh[count] - x_mesh[i]) / h_current) +
                          solve_conju[2 * (i + 1) + 1] * base_func (4, (show_mesh[count] - x_mesh[i]) / h_current);

          show_theta[count] = solve_conju[2 * i] * base_d_func (1, (show_mesh[count] - x_mesh[i]) / h_current) +
                              solve_conju[2 * i + 1] * base_d_func (2, (show_mesh[count] - x_mesh[i]) / h_current) +
                              solve_conju[2 * (i + 1)] * base_d_func (3, (show_mesh[count] - x_mesh[i]) / h_current) +
                              solve_conju[2 * (i + 1) + 1] * base_d_func (4, (show_mesh[count] - x_mesh[i]) / h_current);

          show_moment[count] = data[14] * data[15] *
                              (solve_conju[2 * i] * base_dd_func (1, (show_mesh[count] - x_mesh[i]) / h_current) +
                               solve_conju[2 * i + 1] * base_dd_func (2, (show_mesh[count] - x_mesh[i]) / h_current) +
                               solve_conju[2 * (i + 1)] * base_dd_func (3, (show_mesh[count] - x_mesh[i]) / h_current) +
                               solve_conju[2 * (i + 1) + 1] * base_dd_func (4, (show_mesh[count] - x_mesh[i]) / h_current));   

          show_q[count] = solve_conju[2 * i] * base_ddd_func (1, (show_mesh[count] - x_mesh[i]) / h_current) +
                          solve_conju[2 * i + 1] * base_ddd_func (2, (show_mesh[count] - x_mesh[i]) / h_current) +
                          solve_conju[2 * (i + 1)] * base_ddd_func (3, (show_mesh[count] - x_mesh[i]) / h_current) +
                          solve_conju[2 * (i + 1) + 1] * base_ddd_func (4, (show_mesh[count] - x_mesh[i]) / h_current);                                      

          count++;
        }
    }

/*
  cout << "show_mesh = ";         show_vect (show_mesh, show_count, 1);
  cout << "\n\nshow_w = ";        show_vect (show_w, show_count, 1);
  cout << "\n\nshow_theta = ";    show_vect (show_theta, show_count, 1);
  cout << "\n\nshow_moment = ";   show_vect (show_moment, show_count, 1);
  cout << "\n\nshow_q = ";        show_vect (show_q, show_count, 1);
*/

  write_arrays_to_file("output.txt", show_mesh, show_w, show_theta, show_moment, show_q, show_count);

  system("sage main.sage");

//---------------------------------------------------------------------------------------
/* !!!!!!!!!!!!!!!!! написать нормальную очистку памяти !!!!!!!!!!! */

  delete[] x_mesh;

  return 0;
}